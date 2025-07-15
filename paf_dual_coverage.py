#!/usr/bin/env python3
"""
PAF Dual Coverage Calculator - Query and Target BED Coverage Analysis

This tool analyzes PAF (Pairwise Alignment Format) files to calculate coverage statistics 
for both query and target sequences based on BED-defined regions. It supports multi-core
parallel processing for efficient analysis of large alignment files.

Key Features:
- Dual coverage tracking: separate statistics for query and target sequences
- Multi-core parallel processing with configurable thread count
- Interval tree data structure for efficient region overlap queries
- CIGAR string parsing for accurate base-level coverage
- Progress bar with real-time statistics
- Comprehensive coverage reports with depth distribution

Usage:
    paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed [options]

Author: [Your name]
Version: 1.0
License: MIT
"""

import argparse
import re
import sys
import os
from collections import defaultdict
from typing import List, Tuple, Dict, Set, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import time
import bisect
from dataclasses import dataclass
from tqdm import tqdm


@dataclass
class Interval:
    """
    Represents a genomic interval for interval tree implementation.
    
    This class stores a genomic interval with its start and end positions,
    along with associated data. Used for efficient overlap queries in O(log n + k) time,
    where n is the number of intervals and k is the number of overlapping intervals.
    
    Attributes:
        start: 0-based start position (inclusive)
        end: 0-based end position (exclusive)
        data: Dictionary containing serialized BEDRegion information
    """
    start: int  # 0-based start position
    end: int    # 0-based end position (exclusive)
    data: Dict  # Serialized BEDRegion data for multiprocessing


class IntervalTree:
    """
    Simple interval tree implementation for fast overlap queries.
    
    This implementation uses a sorted list with binary search for efficiency.
    While not a true interval tree (which would use a balanced tree structure),
    it provides good performance for our use case where intervals are loaded
    once and queried many times.
    
    Time Complexity:
    - Add interval: O(1) amortized
    - Query: O(log n + k) where n = total intervals, k = overlapping intervals
    """
    
    def __init__(self):
        """Initialize an empty interval tree."""
        self.intervals = []  # List of Interval objects
        self.sorted = False  # Flag to track if intervals are sorted
    
    def add(self, start: int, end: int, data: Dict) -> None:
        """
        Add an interval to the tree.
        
        The interval is added to an unsorted list for efficiency during loading.
        Sorting is deferred until the first query operation.
        
        Args:
            start: Start position of the interval (0-based, inclusive)
            end: End position of the interval (0-based, exclusive)
            data: Associated data for the interval (must be serializable)
        """
        self.intervals.append(Interval(start, end, data))
        self.sorted = False
    
    def _ensure_sorted(self) -> None:
        """
        Ensure intervals are sorted by start position for binary search.
        
        This method is called before any query operation to maintain
        the sorted invariant required for efficient searching.
        """
        if not self.sorted:
            self.intervals.sort(key=lambda x: x.start)
            self.sorted = True
    
    def query(self, start: int, end: int) -> List[Dict]:
        """
        Find all intervals overlapping with the query range.
        
        Two intervals overlap if:
        - interval.start < query.end AND interval.end > query.start
        
        Algorithm:
        1. Use binary search to find the first interval that could overlap
        2. Iterate through subsequent intervals until no more overlaps possible
        
        Args:
            start: Query start position (0-based, inclusive)
            end: Query end position (0-based, exclusive)
            
        Returns:
            List of data dictionaries for overlapping intervals
        """
        self._ensure_sorted()
        results = []
        
        # Binary search to find the first potentially overlapping interval
        # We search for intervals where interval.end > start (necessary for overlap)
        left = bisect.bisect_left(self.intervals, start, key=lambda x: x.end)
        
        # Check all intervals from this point until we pass the query region
        for i in range(left, len(self.intervals)):
            interval = self.intervals[i]
            
            # If interval starts after query ends, no more overlaps possible
            if interval.start >= end:
                break
                
            # Check for actual overlap
            if interval.end > start:
                results.append(interval.data)
        
        return results


class BEDRegion:
    """
    Represents a genomic region from a BED file.
    
    This class stores information about a genomic region and tracks coverage
    at base-level resolution. It supports both query and target coordinate
    systems for dual coverage analysis.
    
    BED format uses 0-based, half-open intervals:
    - Start is 0-based (first base of chromosome is 0)
    - End is exclusive (points to one base after the last base in the region)
    
    Attributes:
        chrom: Chromosome/sequence name
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        name: Optional region name/identifier
        source: Source BED file identifier (BED1/BED2)
        coord_type: Coordinate type ('query' or 'target')
        length: Region length in bases
        coverage: Dictionary mapping position to coverage depth
        region_id: Unique identifier for this region
    """
    
    def __init__(self, chrom: str, start: int, end: int, name: str = None, 
                 source: str = None, coord_type: str = None):
        """
        Initialize a BED region.
        
        Args:
            chrom: Chromosome/sequence name
            start: Start position (0-based, inclusive)
            end: End position (0-based, exclusive)
            name: Optional region name (defaults to chrom:start-end)
            source: Source BED file identifier (BED1/BED2)
            coord_type: Coordinate type ('query' or 'target')
        """
        self.chrom = chrom
        self.start = start  # 0-based, inclusive
        self.end = end      # 0-based, exclusive
        self.name = name or f"{chrom}:{start}-{end}"
        self.source = source  # BED1 or BED2
        self.coord_type = coord_type  # 'query' or 'target'
        self.length = end - start
        self.coverage = defaultdict(int)  # position -> depth mapping
        
        # Create unique identifier for efficient lookups
        self.region_id = f"{coord_type}:{chrom}:{start}-{end}"
    
    def __str__(self):
        """
        String representation using 1-based coordinates for display.
        
        While BED files use 0-based coordinates internally, we display
        1-based coordinates for human readability.
        """
        return f"{self.chrom}:{self.start+1}-{self.end}"
    
    def to_dict(self) -> Dict:
        """
        Convert to a serializable dictionary for multiprocessing.
        
        The coverage defaultdict cannot be pickled directly, so we only
        include the basic region information. Coverage will be accumulated
        separately during processing.
        
        Returns:
            Dictionary with region information (without coverage data)
        """
        return {
            'chrom': self.chrom,
            'start': self.start,
            'end': self.end,
            'name': self.name,
            'source': self.source,
            'coord_type': self.coord_type,
            'region_id': self.region_id
        }


def parse_bed_file(bed_file: str, source_id: str, coord_type: str) -> Tuple[List[BEDRegion], Dict[str, List]]:
    """
    Parse a BED file and create region objects with interval trees.
    
    This function reads a BED file, creates BEDRegion objects for each entry,
    and prepares interval tree data structures for efficient overlap queries.
    The interval trees are serialized as lists for multiprocessing compatibility.
    
    BED format specification:
    - Column 1: Chromosome name
    - Column 2: Start position (0-based)
    - Column 3: End position (0-based, exclusive)
    - Column 4: Name (optional)
    - Columns 5+: Additional fields (ignored)
    
    Args:
        bed_file: Path to the BED file
        source_id: Source identifier (BED1 for query, BED2 for target)
        coord_type: Coordinate type ('query' or 'target')
        
    Returns:
        Tuple of:
        - List of BEDRegion objects
        - Dictionary mapping chromosome names to serialized interval tree data
        
    Raises:
        SystemExit: If file is not found
    """
    regions = []
    serialized_trees = defaultdict(list)
    
    try:
        with open(bed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                # Split by tab (BED is tab-delimited)
                fields = line.split('\t')
                
                # Validate minimum required fields
                if len(fields) < 3:
                    print(f"Warning: {bed_file} line {line_num} has invalid format, skipping", 
                          file=sys.stderr)
                    continue
                
                # Parse required fields
                chrom = fields[0]
                try:
                    start = int(fields[1])  # 0-based start
                    end = int(fields[2])    # 0-based end (exclusive)
                except ValueError:
                    print(f"Warning: {bed_file} line {line_num} has invalid coordinates, skipping", 
                          file=sys.stderr)
                    continue
                
                # Validate coordinates
                if start < 0 or end < start:
                    print(f"Warning: {bed_file} line {line_num} has invalid interval [{start}, {end}), skipping", 
                          file=sys.stderr)
                    continue
                
                # Parse optional name field (column 4)
                name = fields[3] if len(fields) > 3 else None
                
                # Create region object
                region = BEDRegion(chrom, start, end, name, source_id, coord_type)
                regions.append(region)
                
                # Add to serialized interval tree data
                # Format: (start, end, region_data_dict)
                serialized_trees[chrom].append((start, end, region.to_dict()))
                
    except FileNotFoundError:
        print(f"Error: File not found - {bed_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to parse BED file {bed_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Convert defaultdict to regular dict for pickling
    return regions, dict(serialized_trees)


def parse_cigar(cigar_string: str) -> List[Tuple[int, str]]:
    """
    Parse a CIGAR string into a list of operations.
    
    CIGAR (Compact Idiosyncratic Gapped Alignment Report) format encodes
    the alignment structure as a series of operation codes with lengths.
    
    CIGAR operations and their effects on coverage:
    - M: Alignment match (sequence match or mismatch) - consumes query and target
    - =: Sequence match - consumes query and target
    - X: Sequence mismatch - consumes query and target
    - I: Insertion to target - consumes query only
    - D: Deletion from target - consumes target only
    - N: Skipped region in target (e.g., intron) - consumes target only
    - S: Soft clipping - consumes query only (bases present in sequence)
    - H: Hard clipping - consumes neither (bases not present in sequence)
    - P: Padding - silent deletion from padded reference
    
    Args:
        cigar_string: CIGAR string from PAF file (e.g., "100M2D50M")
        
    Returns:
        List of (length, operation) tuples, e.g., [(100, 'M'), (2, 'D'), (50, 'M')]
    """
    if not cigar_string or cigar_string == '*':
        return []
    
    # Regular expression to parse CIGAR operations
    # Matches: one or more digits followed by a valid operation character
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    matches = pattern.findall(cigar_string)
    
    # Convert string lengths to integers
    return [(int(length), op) for length, op in matches]


def calculate_coverage_from_cigar(query_start: int, target_start: int, 
                                cigar_ops: List[Tuple[int, str]], 
                                strand: str,
                                query_trees: Dict[str, IntervalTree],
                                target_trees: Dict[str, IntervalTree],
                                query_chrom: str,
                                target_chrom: str,
                                query_len: int) -> Tuple[Dict, Dict]:
    """
    Calculate coverage for both query and target sequences based on CIGAR operations.
    
    This function is the core of accurate coverage calculation. It walks through
    the CIGAR string operation by operation, tracking positions in both query
    and target sequences, and updating coverage for any overlapping regions.
    
    Special handling for negative strand alignments:
    - PAF reports query coordinates on the original strand
    - For negative strand, we need to convert to actual genomic coordinates
    - Target coordinates are always on the forward strand
    
    Algorithm:
    1. Initialize positions in query and target
    2. For each CIGAR operation:
       - Determine which sequence(s) are consumed
       - Find overlapping regions using interval trees
       - Update coverage for each overlapping base
       - Advance position counters appropriately
    
    Args:
        query_start: Start position in query (0-based)
        target_start: Start position in target (0-based)
        cigar_ops: List of CIGAR operations from parse_cigar()
        strand: Alignment strand ('+' or '-')
        query_trees: Interval trees for query regions (by chromosome)
        target_trees: Interval trees for target regions (by chromosome)
        query_chrom: Query chromosome/sequence name
        target_chrom: Target chromosome/sequence name
        query_len: Total length of query sequence (needed for negative strand)
        
    Returns:
        Tuple of (query_coverage, target_coverage) where each is a dict:
        - Key: region_id
        - Value: dict mapping position to coverage depth
    """
    # Initialize coverage tracking
    # Using nested defaultdict for efficient updates
    query_coverage = defaultdict(lambda: defaultdict(int))
    target_coverage = defaultdict(lambda: defaultdict(int))
    
    # Current positions in query and target
    query_pos = query_start
    target_pos = target_start
    
    # Get interval trees for the relevant chromosomes
    query_tree = query_trees.get(query_chrom) if query_chrom in query_trees else None
    target_tree = target_trees.get(target_chrom) if target_chrom in target_trees else None
    
    # Process each CIGAR operation in order
    for length, op in cigar_ops:
        if op in ['M', '=', 'X']:  # Match/mismatch operations
            # These operations consume both query and target sequences
            
            # Update query coverage
            if query_tree:
                if strand == '-':
                    # For negative strand alignments, convert coordinates
                    # PAF reports positions on the original strand, but we need
                    # actual genomic coordinates for coverage calculation
                    # 
                    # Example: Query length = 1000, alignment at positions 100-200
                    # On negative strand, this corresponds to positions 800-900
                    actual_query_start = query_len - query_pos - length
                    actual_query_end = query_len - query_pos
                    
                    # Find overlapping regions
                    overlapping_regions = query_tree.query(actual_query_start, actual_query_end)
                    
                    for region_data in overlapping_regions:
                        # Calculate the overlap between alignment and region
                        overlap_start = max(actual_query_start, region_data['start'])
                        overlap_end = min(actual_query_end, region_data['end'])
                        region_id = region_data['region_id']
                        
                        # Update coverage for each base in the overlap
                        for pos in range(overlap_start, overlap_end):
                            query_coverage[region_id][pos] += 1
                else:
                    # Positive strand - use coordinates directly
                    overlapping_regions = query_tree.query(query_pos, query_pos + length)
                    
                    for region_data in overlapping_regions:
                        overlap_start = max(query_pos, region_data['start'])
                        overlap_end = min(query_pos + length, region_data['end'])
                        region_id = region_data['region_id']
                        
                        for pos in range(overlap_start, overlap_end):
                            query_coverage[region_id][pos] += 1
            
            # Update target coverage (strand doesn't affect target coordinates)
            if target_tree:
                overlapping_regions = target_tree.query(target_pos, target_pos + length)
                
                for region_data in overlapping_regions:
                    overlap_start = max(target_pos, region_data['start'])
                    overlap_end = min(target_pos + length, region_data['end'])
                    region_id = region_data['region_id']
                    
                    for pos in range(overlap_start, overlap_end):
                        target_coverage[region_id][pos] += 1
            
            # Advance both positions
            query_pos += length
            target_pos += length
            
        elif op == 'I':  # Insertion to target
            # Insertion means query has extra bases not in target
            # Consumes query but not target
            query_pos += length
            # target_pos stays the same
            
        elif op == 'D':  # Deletion from target
            # Deletion means target has extra bases not in query
            # Consumes target but not query
            # query_pos stays the same
            target_pos += length
            
        elif op == 'N':  # Skipped region in target
            # Used for introns in RNA-to-genome alignments
            # Similar to deletion but implies a large gap
            target_pos += length
            # query_pos stays the same
            
        elif op == 'S':  # Soft clipping
            # Clipped sequences are present in the query but not aligned
            # Consumes query but not target
            query_pos += length
            # target_pos stays the same
            
        elif op == 'H':  # Hard clipping
            # Clipped sequences are not present in the PAF/SAM record
            # Doesn't consume either sequence
            pass
    
    # Convert defaultdicts to regular dicts for pickling
    return dict(query_coverage), dict(target_coverage)


def parse_paf_line(line: str) -> Optional[Dict]:
    """
    Parse a single line from a PAF file.
    
    PAF (Pairwise Alignment Format) is a text format for storing pairwise alignments.
    It's designed to be a compact alternative to SAM for storing alignments.
    
    PAF format specification (tab-delimited):
    Col  Field         Type    Description
    1    query_name    string  Query sequence name
    2    query_len     int     Query sequence length
    3    query_start   int     Query start (0-based, inclusive)
    4    query_end     int     Query end (0-based, exclusive)
    5    strand        char    Relative strand (+ or -)
    6    target_name   string  Target sequence name
    7    target_len    int     Target sequence length
    8    target_start  int     Target start (0-based, inclusive)
    9    target_end    int     Target end (0-based, exclusive)
    10   matches       int     Number of matching bases
    11   block_len     int     Alignment block length
    12   mapq          int     Mapping quality (0-255; 255 for missing)
    13+  tags          string  Optional SAM-style tags
    
    Common optional tags:
    - cg:Z:CIGAR    CIGAR string for detailed alignment
    - NM:i:n        Number of mismatches
    - AS:i:n        Alignment score
    
    Args:
        line: A line from the PAF file
        
    Returns:
        Dictionary with parsed fields, or None if parsing fails
    """
    fields = line.strip().split('\t')
    
    # Validate minimum required fields
    if len(fields) < 12:
        return None
    
    # Parse mandatory fields
    try:
        result = {
            'query_name': fields[0],
            'query_len': int(fields[1]),
            'query_start': int(fields[2]),
            'query_end': int(fields[3]),
            'strand': fields[4],
            'target_name': fields[5],
            'target_len': int(fields[6]),
            'target_start': int(fields[7]),
            'target_end': int(fields[8]),
            'matches': int(fields[9]),
            'alignment_len': int(fields[10]),
            'mapq': int(fields[11])
        }
    except (ValueError, IndexError):
        return None
    
    # Parse optional fields (SAM-style tags)
    for field in fields[12:]:
        # Tags format: TAG:TYPE:VALUE
        if ':' in field:
            parts = field.split(':', 2)
            if len(parts) >= 3:
                tag, tag_type, value = parts[0], parts[1], parts[2]
                
                # Look for CIGAR string
                if tag == 'cg' and tag_type == 'Z':
                    result['cigar'] = value
                # Could parse other tags here if needed
                # elif tag == 'NM' and tag_type == 'i':
                #     result['mismatches'] = int(value)
    
    return result


def process_paf_chunk(chunk_data: Tuple[str, int, int, Dict, Dict, int, int]) -> Dict:
    """
    Process a chunk of the PAF file in a separate process.
    
    This function is designed for parallel execution. Each process handles
    an independent chunk of the PAF file, accumulating coverage statistics
    that will be merged in the main process.
    
    The chunk boundaries are adjusted to complete lines to avoid partial
    record processing.
    
    Algorithm:
    1. Reconstruct interval trees from serialized data
    2. Seek to chunk start position
    3. Process lines until chunk end position
    4. Calculate coverage using CIGAR strings when available
    5. Return accumulated statistics and coverage data
    
    Args:
        chunk_data: Tuple containing:
            - paf_file: Path to PAF file
            - start_pos: Starting byte position in file
            - end_pos: Ending byte position in file  
            - query_trees_data: Serialized query interval trees
            - target_trees_data: Serialized target interval trees
            - min_mapq: Minimum mapping quality threshold
            - chunk_id: Identifier for this chunk (for debugging)
            
    Returns:
        Dictionary containing:
            - stats: Processing statistics for this chunk
            - query_coverage: Query coverage data
            - target_coverage: Target coverage data
            - chunk_id: Chunk identifier
            - query_chroms: List of query chromosomes seen
            - target_chroms: List of target chromosomes seen
    """
    # Unpack arguments
    paf_file, start_pos, end_pos, query_trees_data, target_trees_data, min_mapq, chunk_id = chunk_data
    
    # Rebuild interval trees from serialized data
    # This is necessary because IntervalTree objects can't be pickled directly
    query_trees = {}
    for chrom, intervals in query_trees_data.items():
        tree = IntervalTree()
        for start, end, region_data in intervals:
            tree.add(start, end, region_data)
        query_trees[chrom] = tree
    
    target_trees = {}
    for chrom, intervals in target_trees_data.items():
        tree = IntervalTree()
        for start, end, region_data in intervals:
            tree.add(start, end, region_data)
        target_trees[chrom] = tree
    
    # Initialize local statistics for this chunk
    local_stats = {
        'alignment_count': 0,      # Valid alignments processed
        'filtered_count': 0,       # Alignments filtered by quality
        'no_cigar_count': 0,      # Alignments without CIGAR
        'query_chroms': set(),    # Unique query chromosomes
        'target_chroms': set(),   # Unique target chromosomes  
        'lines_processed': 0      # Total lines read
    }
    
    # Initialize local coverage data
    local_query_coverage = defaultdict(lambda: defaultdict(int))
    local_target_coverage = defaultdict(lambda: defaultdict(int))
    
    try:
        with open(paf_file, 'r') as f:
            # Seek to starting position for this chunk
            f.seek(start_pos)
            
            # Process lines until we reach the end position
            while f.tell() < end_pos:
                line = f.readline()
                if not line:
                    break
                
                local_stats['lines_processed'] += 1
                line = line.strip()
                
                if not line:
                    continue
                
                # Parse PAF record
                paf_record = parse_paf_line(line)
                if not paf_record:
                    continue
                
                # Track chromosomes seen
                local_stats['query_chroms'].add(paf_record['query_name'])
                local_stats['target_chroms'].add(paf_record['target_name'])
                
                # Apply mapping quality filter
                if paf_record['mapq'] < min_mapq:
                    local_stats['filtered_count'] += 1
                    continue
                
                # Calculate coverage based on CIGAR availability
                if 'cigar' not in paf_record:
                    # No CIGAR string - use alignment blocks for coverage
                    # This is less accurate as it doesn't account for indels
                    local_stats['no_cigar_count'] += 1
                    
                    # Simple coverage for query regions
                    if paf_record['query_name'] in query_trees:
                        query_tree = query_trees[paf_record['query_name']]
                        overlapping_regions = query_tree.query(
                            paf_record['query_start'], 
                            paf_record['query_end']
                        )
                        
                        for region_data in overlapping_regions:
                            overlap_start = max(paf_record['query_start'], region_data['start'])
                            overlap_end = min(paf_record['query_end'], region_data['end'])
                            region_id = region_data['region_id']
                            
                            # Update coverage for the overlapping region
                            for pos in range(overlap_start, overlap_end):
                                local_query_coverage[region_id][pos] += 1
                    
                    # Simple coverage for target regions
                    if paf_record['target_name'] in target_trees:
                        target_tree = target_trees[paf_record['target_name']]
                        overlapping_regions = target_tree.query(
                            paf_record['target_start'], 
                            paf_record['target_end']
                        )
                        
                        for region_data in overlapping_regions:
                            overlap_start = max(paf_record['target_start'], region_data['start'])
                            overlap_end = min(paf_record['target_end'], region_data['end'])
                            region_id = region_data['region_id']
                            
                            for pos in range(overlap_start, overlap_end):
                                local_target_coverage[region_id][pos] += 1
                else:
                    # CIGAR-based coverage calculation (more accurate)
                    cigar_ops = parse_cigar(paf_record['cigar'])
                    query_cov, target_cov = calculate_coverage_from_cigar(
                        paf_record['query_start'],
                        paf_record['target_start'],
                        cigar_ops,
                        paf_record['strand'],
                        query_trees,
                        target_trees,
                        paf_record['query_name'],
                        paf_record['target_name'],
                        paf_record['query_len']
                    )
                    
                    # Merge coverage data
                    for region_id, positions in query_cov.items():
                        for pos, depth in positions.items():
                            local_query_coverage[region_id][pos] += depth
                    
                    for region_id, positions in target_cov.items():
                        for pos, depth in positions.items():
                            local_target_coverage[region_id][pos] += depth
                
                local_stats['alignment_count'] += 1
                
    except Exception as e:
        print(f"Chunk {chunk_id} processing error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    
    # Return results as regular dictionaries (not defaultdicts) for pickling
    return {
        'stats': local_stats,
        'query_coverage': {k: dict(v) for k, v in local_query_coverage.items()},
        'target_coverage': {k: dict(v) for k, v in local_target_coverage.items()},
        'chunk_id': chunk_id,
        'query_chroms': list(local_stats['query_chroms']),
        'target_chroms': list(local_stats['target_chroms'])
    }


def get_file_chunks(file_path: str, num_chunks: int) -> List[Tuple[int, int]]:
    """
    Divide a file into chunks for parallel processing.
    
    This function calculates byte positions that divide the file into
    approximately equal chunks while respecting line boundaries. This
    ensures that no lines are split between chunks.
    
    Algorithm:
    1. Calculate ideal chunk size based on file size
    2. For each chunk, adjust start position to beginning of next line
    3. Adjust end position to end of current line
    4. Handle edge cases (empty file, small files)
    
    Args:
        file_path: Path to the file to be chunked
        num_chunks: Desired number of chunks
        
    Returns:
        List of (start_pos, end_pos) tuples defining each chunk
    """
    file_size = os.path.getsize(file_path)
    
    # Handle empty file
    if file_size == 0:
        return [(0, 0)]
    
    # Calculate ideal chunk size
    chunk_size = max(1, file_size // num_chunks)
    
    chunks = []
    with open(file_path, 'rb') as f:
        for i in range(num_chunks):
            # Calculate nominal start position
            start = i * chunk_size
            
            # Adjust start position to beginning of a line
            # (except for the first chunk)
            if i > 0:
                f.seek(start)
                f.readline()  # Skip potentially incomplete line
                start = f.tell()
            
            # Calculate end position
            if i == num_chunks - 1:
                # Last chunk goes to end of file
                end = file_size
            else:
                # Find end of line at nominal chunk boundary
                end = (i + 1) * chunk_size
                f.seek(end)
                f.readline()  # Move to end of current line
                end = f.tell()
            
            # Only add non-empty chunks
            if start < file_size:
                chunks.append((start, min(end, file_size)))
    
    return chunks


def process_paf_file_parallel(paf_file: str, 
                             query_regions: List[BEDRegion],
                             target_regions: List[BEDRegion],
                             query_trees: Dict[str, List],
                             target_trees: Dict[str, List],
                             min_mapq: int = 0, 
                             num_processes: int = None) -> Dict:
    """
    Process a PAF file using parallel processing.
    
    This function coordinates multiple processes to analyze a large PAF file
    efficiently. It divides the file into chunks, processes them in parallel,
    and merges the results.
    
    The parallel processing strategy:
    1. Divide file into N chunks (one per process)
    2. Each process independently calculates coverage for its chunk
    3. Main process merges results from all chunks
    4. Coverage data is accumulated into the original region objects
    
    Memory considerations:
    - Interval trees are duplicated in each process
    - Coverage data is accumulated in memory
    - For very large files, consider processing in batches
    
    Args:
        paf_file: Path to PAF file
        query_regions: List of query BED regions
        target_regions: List of target BED regions
        query_trees: Serialized query interval trees
        target_trees: Serialized target interval trees
        min_mapq: Minimum mapping quality threshold
        num_processes: Number of processes (None = CPU count)
        
    Returns:
        Dictionary with global statistics
    """
    # Determine number of processes
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()
    
    print(f"Using {num_processes} processes for analysis...", file=sys.stderr)
    
    # Create mapping from region_id to region object for result merging
    region_map = {}
    for region in query_regions + target_regions:
        region_map[region.region_id] = region
    
    # Divide file into chunks
    chunks = get_file_chunks(paf_file, num_processes)
    if not chunks:
        return {
            'alignment_count': 0,
            'filtered_count': 0,
            'no_cigar_count': 0,
            'query_chroms': set(),
            'target_chroms': set()
        }
    
    # Prepare task data for each process
    tasks = [
        (paf_file, start, end, query_trees, target_trees, min_mapq, i)
        for i, (start, end) in enumerate(chunks)
    ]
    
    # Initialize global statistics
    global_stats = {
        'alignment_count': 0,
        'filtered_count': 0,
        'no_cigar_count': 0,
        'query_chroms': set(),
        'target_chroms': set()
    }
    
    # Estimate total lines for progress bar
    # Assuming average line length of 200 bytes (conservative estimate)
    file_size = os.path.getsize(paf_file)
    estimated_lines = file_size // 200
    
    # Process chunks in parallel
    start_time = time.time()
    
    print("Starting process pool...", file=sys.stderr)
    
    # Create progress bar
    with tqdm(total=estimated_lines, desc="Processing PAF file", unit="lines") as pbar:
        # Create process pool
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            # Submit all tasks to the process pool
            future_to_chunk = {executor.submit(process_paf_chunk, task): task[6] 
                              for task in tasks}
            
            # Process results as they complete
            for future in as_completed(future_to_chunk):
                chunk_id = future_to_chunk[future]
                try:
                    result = future.result()
                    
                    # Merge statistics
                    global_stats['alignment_count'] += result['stats']['alignment_count']
                    global_stats['filtered_count'] += result['stats']['filtered_count']
                    global_stats['no_cigar_count'] += result['stats']['no_cigar_count']
                    global_stats['query_chroms'].update(result['query_chroms'])
                    global_stats['target_chroms'].update(result['target_chroms'])
                    
                    # Update progress bar
                    lines_in_chunk = result['stats']['lines_processed']
                    pbar.update(lines_in_chunk)
                    pbar.set_postfix({
                        'alignments': f"{global_stats['alignment_count']:,}",
                        'process': f"{chunk_id+1}/{num_processes}"
                    })
                    
                    # Merge coverage information into region objects
                    # This is the critical step where parallel results are combined
                    
                    # Merge query coverage
                    for region_id, positions in result['query_coverage'].items():
                        if region_id in region_map:
                            region = region_map[region_id]
                            for pos, depth in positions.items():
                                region.coverage[pos] += depth
                    
                    # Merge target coverage
                    for region_id, positions in result['target_coverage'].items():
                        if region_id in region_map:
                            region = region_map[region_id]
                            for pos, depth in positions.items():
                                region.coverage[pos] += depth
                                
                except Exception as e:
                    print(f"\nChunk {chunk_id} processing failed: {e}", file=sys.stderr)
                    import traceback
                    traceback.print_exc()
    
    # Calculate and display performance statistics
    elapsed = time.time() - start_time
    print(f"\nProcessing complete. Time elapsed: {elapsed:.2f} seconds", file=sys.stderr)
    if global_stats['alignment_count'] > 0:
        print(f"Processing speed: {global_stats['alignment_count'] / elapsed:.0f} alignments/second", 
              file=sys.stderr)
    
    return global_stats


def calculate_region_stats(region: BEDRegion) -> Dict:
    """
    Calculate comprehensive coverage statistics for a single region.
    
    This function analyzes the base-level coverage data to compute various
    statistics that help assess the quality and completeness of coverage.
    
    Computed statistics:
    - Coverage percentage: Fraction of bases with at least 1x coverage
    - Average depth: Mean coverage across entire region
    - Average depth (covered only): Mean coverage excluding zero-coverage bases
    - Maximum depth: Highest coverage observed
    - Depth distribution: Histogram of coverage depths
    
    Args:
        region: BEDRegion object with populated coverage data
        
    Returns:
        Dictionary containing:
        - region: The BEDRegion object
        - covered_bases: Number of bases with coverage > 0
        - uncovered_bases: Number of bases with coverage = 0
        - coverage_percent: Percentage of region covered
        - average_depth: Mean coverage across entire region
        - average_depth_covered_only: Mean coverage of covered bases only
        - max_depth: Maximum coverage depth observed
        - depth_distribution: Dict mapping depth to base count
    """
    # Initialize counters
    covered_bases = 0      # Bases with coverage > 0
    total_coverage = 0     # Sum of all coverage values
    max_depth = 0          # Maximum coverage observed
    depth_distribution = defaultdict(int)  # Histogram of depths
    
    # Analyze each position in the region
    for pos in range(region.start, region.end):
        depth = region.coverage.get(pos, 0)
        
        if depth > 0:
            covered_bases += 1
            total_coverage += depth
            max_depth = max(max_depth, depth)
        
        # Build depth distribution
        depth_distribution[depth] += 1
    
    # Calculate summary statistics
    coverage_percent = (covered_bases / region.length * 100) if region.length > 0 else 0
    avg_depth = total_coverage / region.length if region.length > 0 else 0
    avg_depth_covered = total_coverage / covered_bases if covered_bases > 0 else 0
    
    return {
        'region': region,
        'covered_bases': covered_bases,
        'uncovered_bases': region.length - covered_bases,
        'coverage_percent': coverage_percent,
        'average_depth': avg_depth,
        'average_depth_covered_only': avg_depth_covered,
        'max_depth': max_depth,
        'depth_distribution': dict(depth_distribution)
    }


def print_summary(title: str, stats_list: List[Dict], process_stats: Dict, 
                 coord_type: str = None, output_file=sys.stdout):
    """
    Print summary statistics for a set of regions.
    
    This function generates a comprehensive summary report including:
    - Processing statistics (alignments, filtering)
    - Overall coverage metrics
    - Coverage distribution across regions
    - Depth statistics
    
    Args:
        title: Title for the summary section
        stats_list: List of region statistics from calculate_region_stats()
        process_stats: Global processing statistics
        coord_type: Type of coordinates ('query' or 'target')
        output_file: File handle for output (default: stdout)
    """
    print(f"\n=== {title} ===", file=output_file)
    
    # Basic statistics
    print(f"Number of regions: {len(stats_list)}", file=output_file)
    print(f"Alignments processed: {process_stats['alignment_count']:,}", file=output_file)
    print(f"Low quality filtered: {process_stats['filtered_count']:,}", file=output_file)
    print(f"No CIGAR information: {process_stats['no_cigar_count']:,}", file=output_file)
    
    # Show relevant chromosome/sequence information
    if coord_type == 'query':
        chroms = process_stats.get('query_chroms', set())
        print(f"Query sequences involved: {len(chroms)}", file=output_file)
    elif coord_type == 'target':
        chroms = process_stats.get('target_chroms', set())
        print(f"Target sequences involved: {len(chroms)}", file=output_file)
    
    # Calculate overall statistics
    total_length = sum(s['region'].length for s in stats_list)
    total_covered = sum(s['covered_bases'] for s in stats_list)
    fully_covered = sum(1 for s in stats_list if s['coverage_percent'] == 100)
    partially_covered = sum(1 for s in stats_list if 0 < s['coverage_percent'] < 100)
    not_covered = sum(1 for s in stats_list if s['coverage_percent'] == 0)
    
    # Coverage summary
    print(f"\nTotal length: {total_length:,} bp", file=output_file)
    print(f"Covered bases: {total_covered:,} bp", file=output_file)
    if total_length > 0:
        print(f"Overall coverage: {total_covered/total_length*100:.2f}%", file=output_file)
    
    # Region coverage breakdown
    print(f"\nRegion coverage breakdown:", file=output_file)
    print(f"  Fully covered (100%): {fully_covered} regions", file=output_file)
    print(f"  Partially covered (0-100%): {partially_covered} regions", file=output_file)
    print(f"  Not covered (0%): {not_covered} regions", file=output_file)
    
    # Depth statistics
    avg_depths = [s['average_depth'] for s in stats_list if s['average_depth'] > 0]
    if avg_depths:
        print(f"\nAverage depth statistics:", file=output_file)
        print(f"  Mean: {sum(avg_depths)/len(avg_depths):.2f}×", file=output_file)
        print(f"  Min: {min(avg_depths):.2f}×", file=output_file)
        print(f"  Max: {max(avg_depths):.2f}×", file=output_file)


def print_region_stats(stats: Dict, coord_type: str = "", output_file=sys.stdout):
    """
    Print detailed statistics for a single region.
    
    Provides a detailed view of coverage for an individual region,
    useful for identifying specific problem areas or validating coverage.
    
    Args:
        stats: Region statistics from calculate_region_stats()
        coord_type: Coordinate type label ('Query' or 'Target')
        output_file: File handle for output
    """
    region = stats['region']
    
    print(f"\n[{coord_type}] {region.name}", file=output_file)
    print(f"  Location: {region}", file=output_file)
    print(f"  Length: {region.length:,} bp", file=output_file)
    print(f"  Coverage: {stats['covered_bases']:,} bp ({stats['coverage_percent']:.2f}%)", file=output_file)
    print(f"  Average depth: {stats['average_depth']:.2f}×", file=output_file)
    if stats['covered_bases'] > 0:
        print(f"  Average depth (covered only): {stats['average_depth_covered_only']:.2f}×", file=output_file)
    print(f"  Maximum depth: {stats['max_depth']}×", file=output_file)


def main():
    """
    Main function - parse command line arguments and coordinate analysis.
    
    This function:
    1. Parses command line arguments
    2. Validates input files
    3. Loads BED regions and builds interval trees
    4. Processes PAF file in parallel
    5. Calculates coverage statistics
    6. Generates output reports
    """
    # Create argument parser with detailed help
    parser = argparse.ArgumentParser(
        description='Analyze PAF alignment coverage for query and target BED regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Description:
  This tool analyzes PAF alignment files to calculate coverage statistics
  for both query and target sequences based on BED-defined regions.
  
  It uses parallel processing to efficiently handle large alignment files
  and provides comprehensive coverage reports including depth distribution.
  
  BED1: Regions on query sequences
  BED2: Regions on target (reference) sequences

Examples:
  # Basic usage
  %(prog)s alignment.paf query_regions.bed target_regions.bed
  
  # Use 40 processes for large files
  %(prog)s alignment.paf query_regions.bed target_regions.bed -p 40
  
  # Filter low quality alignments and output detailed information
  %(prog)s alignment.paf query_regions.bed target_regions.bed \\
      --min-mapq 30 --output-per-region --output-file results.txt
  
  # Process from stdin (e.g., from minimap2 pipeline)
  minimap2 -x asm5 ref.fa query.fa | %(prog)s - query.bed target.bed

Notes:
  - BED files use 0-based, half-open intervals
  - PAF coordinates are 0-based
  - CIGAR strings provide more accurate coverage than alignment blocks
  - Negative strand alignments are handled correctly
        """
    )
    
    # Positional arguments
    parser.add_argument('paf_file', 
                       help='PAF alignment file (use "-" for stdin)')
    parser.add_argument('query_bed', 
                       help='BED file for query sequences (BED1)')
    parser.add_argument('target_bed', 
                       help='BED file for target sequences (BED2)')
    
    # Optional arguments
    parser.add_argument('--min-mapq', type=int, default=0, 
                       help='Minimum mapping quality (default: 0)')
    parser.add_argument('--processes', '-p', type=int, default=None,
                       help='Number of processes (default: CPU count)')
    parser.add_argument('--output-per-region', action='store_true',
                       help='Output detailed statistics for each region')
    parser.add_argument('--output-uncovered', action='store_true',
                       help='Output detailed information about uncovered regions')
    parser.add_argument('--output-file', '-o', type=str, default=None,
                       help='Output file path (default: stdout)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate input files
    for f in [args.paf_file, args.query_bed, args.target_bed]:
        if f != '-' and not os.path.exists(f):
            print(f"Error: File does not exist - {f}", file=sys.stderr)
            sys.exit(1)
    
    # Handle stdin for PAF file
    if args.paf_file == '-':
        # Create temporary file for stdin content
        import tempfile
        temp_paf = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.paf')
        print("Reading PAF data from stdin...", file=sys.stderr)
        for line in sys.stdin:
            temp_paf.write(line)
        temp_paf.close()
        args.paf_file = temp_paf.name
    
    # Set up output file
    output_file = sys.stdout
    if args.output_file:
        try:
            output_file = open(args.output_file, 'w')
        except IOError as e:
            print(f"Error: Cannot open output file - {e}", file=sys.stderr)
            sys.exit(1)
    
    try:
        # Parse BED files
        print("Reading BED files...", file=sys.stderr)
        
        # Query BED (BED1)
        query_regions, query_trees = parse_bed_file(args.query_bed, 'BED1', 'query')
        print(f"Query BED (BED1): {len(query_regions)} regions", file=sys.stderr)
        
        # Target BED (BED2)  
        target_regions, target_trees = parse_bed_file(args.target_bed, 'BED2', 'target')
        print(f"Target BED (BED2): {len(target_regions)} regions", file=sys.stderr)
        
        # Display system information
        cpu_count = multiprocessing.cpu_count()
        print(f"\nSystem CPU cores: {cpu_count}", file=sys.stderr)
        if args.processes and args.processes > cpu_count:
            print(f"Warning: Requested processes ({args.processes}) exceeds CPU cores ({cpu_count})", 
                  file=sys.stderr)
        
        # Process PAF file
        print("\nProcessing PAF file...", file=sys.stderr)
        process_stats = process_paf_file_parallel(
            args.paf_file, 
            query_regions,
            target_regions,
            query_trees,
            target_trees,
            args.min_mapq, 
            args.processes
        )
        
        # Calculate coverage statistics
        print("\nCalculating coverage statistics...", file=sys.stderr)
        
        # Query statistics
        query_stats = [calculate_region_stats(region) for region in query_regions]
        print_summary("Query Coverage Statistics (BED1)", query_stats, process_stats, 
                     coord_type='query', output_file=output_file)
        
        # Target statistics
        target_stats = [calculate_region_stats(region) for region in target_regions]
        print_summary("Target Coverage Statistics (BED2)", target_stats, process_stats,
                     coord_type='target', output_file=output_file)
        
        # Output per-region statistics if requested
        if args.output_per_region:
            print("\n=== Query Region Details (BED1) ===", file=output_file)
            for stats in query_stats:
                print_region_stats(stats, 'Query', output_file=output_file)
            
            print("\n=== Target Region Details (BED2) ===", file=output_file)
            for stats in target_stats:
                print_region_stats(stats, 'Target', output_file=output_file)
        
        # Output uncovered regions if requested
        if args.output_uncovered:
            print("\n=== Uncovered Regions ===", file=output_file)
            
            # Query uncovered regions
            query_uncovered = [s['region'] for s in query_stats if s['coverage_percent'] == 0]
            if query_uncovered:
                print(f"\nQuery uncovered regions (BED1, {len(query_uncovered)} regions):", 
                      file=output_file)
                for region in query_uncovered:
                    print(f"  {region.name}: {region}", file=output_file)
            
            # Target uncovered regions
            target_uncovered = [s['region'] for s in target_stats if s['coverage_percent'] == 0]
            if target_uncovered:
                print(f"\nTarget uncovered regions (BED2, {len(target_uncovered)} regions):", 
                      file=output_file)
                for region in target_uncovered:
                    print(f"  {region.name}: {region}", file=output_file)
            
            if not query_uncovered and not target_uncovered:
                print("All regions have coverage", file=output_file)
                
    finally:
        # Clean up
        if args.output_file:
            output_file.close()
        
        # Remove temporary file if we created one for stdin
        if args.paf_file.endswith('.paf') and 'temp_paf' in locals():
            os.unlink(args.paf_file)


if __name__ == '__main__':
    # Set multiprocessing start method for compatibility
    # 'spawn' is more portable across platforms than 'fork'
    multiprocessing.set_start_method('spawn', force=True)
    main()
