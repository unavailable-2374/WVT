#!/usr/bin/env python3
"""
PAF Coverage Calculator with BED Support

A general-purpose tool for calculating base-level coverage statistics from PAF
(Pairwise Alignment Format) files against BED-defined regions. This tool can 
analyze coverage for any reference sequences and supports merging multiple BED files.

Key Features:
- Base-level coverage calculation using CIGAR strings
- Support for multiple BED files with optional merging
- Mapping quality filtering
- Detailed per-region statistics
- Identification of uncovered regions
- Simple single-threaded implementation

Usage:
    paf_coverage.py alignment.paf regions1.bed regions2.bed [options]

Author: [Your name]
Version: 1.0
License: MIT
"""

import argparse
import re
import sys
from collections import defaultdict
from typing import List, Tuple, Dict, Set, Optional


class BEDRegion:
    """
    Represents a genomic region from a BED file.
    
    This class encapsulates a genomic region with its coordinates and tracks
    coverage information at base-level resolution. The coverage data allows
    for accurate calculation of coverage depth and breadth statistics.
    
    BED coordinate system:
    - 0-based: First base of chromosome is position 0
    - Half-open intervals: [start, end) where end is exclusive
    - Example: First 100 bases = [0, 100)
    
    Attributes:
        chrom: Chromosome or sequence name
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        name: Region name/identifier
        length: Region length in bases
        coverage: Dictionary mapping position to coverage depth
    """
    
    def __init__(self, chrom: str, start: int, end: int, name: str = None):
        """
        Initialize a BED region.
        
        Args:
            chrom: Chromosome or sequence name
            start: Start position (0-based, inclusive)
            end: End position (0-based, exclusive)
            name: Optional region name/identifier
        """
        self.chrom = chrom
        self.start = start  # 0-based, inclusive
        self.end = end      # 0-based, exclusive
        self.name = name or f"{chrom}:{start}-{end}"
        self.length = end - start
        
        # Dictionary to track coverage depth at each position
        # Key: genomic position (0-based)
        # Value: coverage depth (number of overlapping alignments)
        self.coverage = defaultdict(int)
    
    def __str__(self):
        """
        Return string representation with 1-based coordinates for display.
        
        While internal representation uses 0-based coordinates (following BED
        convention), we display 1-based coordinates for human readability.
        """
        return f"{self.chrom}:{self.start+1}-{self.end}"


def parse_bed_file(bed_file: str) -> List[BEDRegion]:
    """
    Parse a BED file and create BEDRegion objects.
    
    This function reads a standard BED format file and creates BEDRegion
    objects for each valid entry. It handles common BED file issues like
    comments, empty lines, and invalid formats.
    
    BED format specification (tab-delimited):
    Column 1: chrom - Chromosome name
    Column 2: chromStart - Start position (0-based)
    Column 3: chromEnd - End position (0-based, exclusive)
    Column 4: name - Feature name (optional)
    Column 5+: Additional fields (ignored)
    
    Example BED file:
    chr1    1000    2000    gene1
    chr1    5000    6000    gene2
    chr2    1000    3000    gene3
    
    Args:
        bed_file: Path to the BED file
        
    Returns:
        List of BEDRegion objects
        
    Raises:
        SystemExit: If file is not found or has critical errors
    """
    regions = []
    
    try:
        with open(bed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Skip empty lines and comments
                # BED files often have header lines starting with #
                if not line or line.startswith('#'):
                    continue
                
                # Split by tab (BED is tab-delimited)
                fields = line.split('\t')
                
                # Validate minimum required fields
                if len(fields) < 3:
                    print(f"Warning: {bed_file} line {line_num} has invalid format, skipping", 
                          file=sys.stderr)
                    continue
                
                # Parse chromosome name
                chrom = fields[0]
                
                # Parse coordinates
                try:
                    start = int(fields[1])  # 0-based start
                    end = int(fields[2])    # 0-based end (exclusive)
                except ValueError:
                    print(f"Warning: {bed_file} line {line_num} has invalid coordinates, skipping", 
                          file=sys.stderr)
                    continue
                
                # Validate coordinates
                if start < 0:
                    print(f"Warning: {bed_file} line {line_num} has negative start position, skipping", 
                          file=sys.stderr)
                    continue
                    
                if end <= start:
                    print(f"Warning: {bed_file} line {line_num} has invalid interval (end <= start), skipping", 
                          file=sys.stderr)
                    continue
                
                # Parse optional name field (4th column)
                name = fields[3] if len(fields) > 3 else None
                
                # Create and store region
                regions.append(BEDRegion(chrom, start, end, name))
                
    except FileNotFoundError:
        print(f"Error: File not found - {bed_file}", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error: Failed to read file {bed_file} - {e}", file=sys.stderr)
        sys.exit(1)
    
    # Report statistics
    print(f"Loaded {len(regions)} regions from {bed_file}", file=sys.stderr)
    
    return regions


def parse_cigar(cigar_string: str) -> List[Tuple[int, str]]:
    """
    Parse a CIGAR string into a list of operations.
    
    CIGAR (Compact Idiosyncratic Gapped Alignment Report) format represents
    how a sequence aligns to a reference using a compact string notation.
    Each operation is encoded as a length followed by an operation code.
    
    Example: "100M2D50M1I25M"
    Means: 100 matches, 2 deletions, 50 matches, 1 insertion, 25 matches
    
    CIGAR operations and their effects:
    - M: Alignment match (can be sequence match or mismatch)
         Consumes: query and reference
    - =: Sequence match (exact match)
         Consumes: query and reference
    - X: Sequence mismatch
         Consumes: query and reference
    - I: Insertion to the reference
         Consumes: query only
    - D: Deletion from the reference
         Consumes: reference only
    - N: Skipped region from the reference (e.g., intron)
         Consumes: reference only
    - S: Soft clipping (clipped sequences present in SEQ)
         Consumes: query only
    - H: Hard clipping (clipped sequences NOT present in SEQ)
         Consumes: neither
    - P: Padding (silent deletion from padded reference)
         Consumes: neither
    
    Args:
        cigar_string: CIGAR string from PAF file
        
    Returns:
        List of (length, operation) tuples
        Example: [(100, 'M'), (2, 'D'), (50, 'M')]
    """
    # Handle missing CIGAR
    if not cigar_string or cigar_string == '*':
        return []
    
    # Regular expression to match CIGAR operations
    # Pattern: one or more digits followed by an operation character
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    
    # Find all matches in the CIGAR string
    matches = pattern.findall(cigar_string)
    
    # Convert to list of tuples with integer lengths
    return [(int(length), op) for length, op in matches]


def calculate_coverage_from_cigar(ref_start: int, cigar_ops: List[Tuple[int, str]], 
                                regions: List[BEDRegion], chrom: str) -> None:
    """
    Calculate coverage based on CIGAR operations.
    
    This function is the core of accurate coverage calculation. It walks through
    the CIGAR string operation by operation, tracking the reference position
    and updating coverage for any regions that overlap with aligned segments.
    
    The algorithm:
    1. Start at the given reference position
    2. For each CIGAR operation:
       - If it consumes reference bases (M, =, X, D, N), check for region overlaps
       - Update coverage for overlapping bases
       - Advance the reference position accordingly
    3. Operations that don't consume reference (I, S, H) don't affect coverage
    
    Why CIGAR-based coverage is important:
    - Simple alignment blocks don't account for insertions/deletions
    - CIGAR provides base-level accuracy
    - Essential for detecting partial coverage of regions
    
    Args:
        ref_start: Starting position on reference (0-based)
        cigar_ops: List of CIGAR operations from parse_cigar()
        regions: List of BED regions to check for overlap
        chrom: Chromosome/sequence name
        
    Note:
        This function modifies the coverage dictionaries of overlapping regions
        in-place for efficiency.
    """
    # Current position on the reference sequence
    ref_pos = ref_start
    
    # Process each CIGAR operation sequentially
    for length, op in cigar_ops:
        
        if op in ['M', '=', 'X']:  # Match/mismatch operations
            # These operations represent aligned bases that consume both
            # query and reference sequences. They contribute to coverage.
            
            # Find all regions that overlap with this alignment segment
            for region in regions:
                # Skip regions on different chromosomes
                if region.chrom != chrom:
                    continue
                
                # Calculate the overlap between alignment segment and region
                # Overlap exists if: segment_start < region_end AND segment_end > region_start
                overlap_start = max(ref_pos, region.start)
                overlap_end = min(ref_pos + length, region.end)
                
                if overlap_start < overlap_end:
                    # We have an overlap! Update coverage for each base
                    for pos in range(overlap_start, overlap_end):
                        region.coverage[pos] += 1
            
            # Advance reference position
            ref_pos += length
            
        elif op == 'D':  # Deletion from reference
            # Deletion means the reference has bases that are not in the query
            # These bases are "covered" by the alignment gap
            # Whether to count deletions as coverage depends on use case
            # Here we skip them (no coverage update) but advance position
            ref_pos += length
            
        elif op == 'N':  # Skipped region (e.g., intron in RNA-seq)
            # Similar to deletion but represents expected gaps (like introns)
            # Typically not counted as coverage
            ref_pos += length
            
        elif op == 'I':  # Insertion to reference
            # Insertion means the query has extra bases not in reference
            # Doesn't consume reference positions, so no position advance
            # Reference position stays the same
            pass
            
        elif op in ['S', 'H']:  # Soft/hard clipping
            # Clipping represents unaligned query sequences
            # S (soft): sequences are present in the record but not aligned
            # H (hard): sequences are not present in the record
            # Neither affects reference position or coverage
            pass
        
        elif op == 'P':  # Padding
            # Padding is used in multiple sequence alignments
            # Represents a silent deletion from padded reference
            # Rarely seen in pairwise alignments
            pass


def parse_paf_line(line: str) -> Optional[Dict]:
    """
    Parse a single line from a PAF file.
    
    PAF (Pairwise Alignment Format) is a text format describing sequence alignments.
    It was designed as a compact alternative to SAM format, especially for
    long-read alignments and genome-to-genome alignments.
    
    PAF format specification (tab-delimited, 12+ columns):
    
    Mandatory fields (columns 1-12):
    Col  Field           Type    Description
    1    query_name      string  Query sequence name
    2    query_len       int     Query sequence total length
    3    query_start     int     Query start (0-based, inclusive)
    4    query_end       int     Query end (0-based, exclusive)
    5    strand          char    Relative strand: '+' (forward) or '-' (reverse)
    6    target_name     string  Target sequence name
    7    target_len      int     Target sequence total length
    8    target_start    int     Target start on original strand (0-based)
    9    target_end      int     Target end on original strand (0-based, exclusive)
    10   num_matches     int     Number of matching bases in the alignment
    11   alignment_len   int     Alignment block length (including gaps)
    12   mapping_qual    int     Mapping quality (0-255; 255 for missing)
    
    Optional fields (columns 13+):
    SAM-style typed key-value pairs in TAG:TYPE:VALUE format
    
    Common tags:
    - cg:Z:CIGAR    CIGAR string of the alignment
    - NM:i:n        Total number of mismatches and gaps
    - AS:i:n        DP alignment score
    - tp:A:P/S      Type of alignment (P=primary, S=secondary)
    - cm:i:n        Number of minimizers on the chain
    - s1:i:n        Chaining score
    - dv:f:n        Divergence (1 - identity)
    
    Args:
        line: A line from the PAF file
        
    Returns:
        Dictionary with parsed fields, or None if parsing fails
    """
    # Split line by tabs
    fields = line.strip().split('\t')
    
    # Validate minimum required fields
    if len(fields) < 12:
        return None
    
    # Parse mandatory fields with error handling
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
        # Invalid format or non-numeric values
        return None
    
    # Validate parsed values
    if result['query_start'] < 0 or result['query_end'] < result['query_start']:
        return None
    if result['target_start'] < 0 or result['target_end'] < result['target_start']:
        return None
    if result['strand'] not in ['+', '-']:
        return None
    if result['mapq'] < 0 or result['mapq'] > 255:
        return None
    
    # Parse optional fields (SAM-style tags)
    for field in fields[12:]:
        # Tags should have format TAG:TYPE:VALUE
        if ':' not in field:
            continue
            
        parts = field.split(':', 2)
        if len(parts) < 3:
            continue
            
        tag, tag_type, value = parts
        
        # Parse specific tags we're interested in
        if tag == 'cg' and tag_type == 'Z':
            # CIGAR string
            result['cigar'] = value
        elif tag == 'NM' and tag_type == 'i':
            # Number of mismatches/gaps
            try:
                result['mismatches'] = int(value)
            except ValueError:
                pass
        elif tag == 'AS' and tag_type == 'i':
            # Alignment score
            try:
                result['align_score'] = int(value)
            except ValueError:
                pass
        # Additional tags can be parsed here as needed
    
    return result


def calculate_region_stats(region: BEDRegion) -> Dict:
    """
    Calculate comprehensive coverage statistics for a single region.
    
    This function analyzes the base-level coverage data collected during
    PAF processing to generate various statistics that characterize the
    coverage quality and completeness.
    
    Calculated metrics:
    1. Coverage breadth: Percentage of bases with at least 1x coverage
    2. Coverage depth: Average number of reads covering each base
    3. Coverage uniformity: Distribution of coverage depths
    4. Coverage gaps: Identification of zero-coverage areas
    
    The statistics help answer questions like:
    - Is the entire region covered?
    - What is the average sequencing depth?
    - Are there coverage gaps or low-coverage areas?
    - Is coverage uniform or highly variable?
    
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
        - depth_distribution: Histogram of coverage depths
    """
    # Initialize counters and statistics
    covered_bases = 0        # Bases with at least 1x coverage
    total_coverage = 0       # Sum of all coverage values
    max_depth = 0           # Maximum depth observed
    depth_distribution = defaultdict(int)  # Histogram: depth -> count
    
    # Analyze each position in the region
    for pos in range(region.start, region.end):
        # Get coverage depth at this position (0 if not in coverage dict)
        depth = region.coverage.get(pos, 0)
        
        # Update statistics
        if depth > 0:
            covered_bases += 1
            total_coverage += depth
            max_depth = max(max_depth, depth)
        
        # Build depth distribution for histogram
        depth_distribution[depth] += 1
    
    # Calculate derived statistics
    # Coverage breadth: what fraction of the region is covered?
    coverage_percent = (covered_bases / region.length * 100) if region.length > 0 else 0
    
    # Average depth across the entire region (including zero-coverage bases)
    avg_depth = total_coverage / region.length if region.length > 0 else 0
    
    # Average depth considering only covered bases
    # This shows the typical depth where we do have coverage
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


def process_paf_file(paf_file: str, regions: List[BEDRegion], min_mapq: int = 0) -> Dict:
    """
    Process a PAF file and calculate coverage for the given regions.
    
    This is the main processing function that reads through a PAF file,
    filters alignments based on quality, and updates coverage information
    for all overlapping regions.
    
    Processing steps:
    1. Organize regions by chromosome for efficient lookup
    2. Read PAF file line by line
    3. Parse each alignment and apply quality filters
    4. Calculate coverage using CIGAR strings when available
    5. Track processing statistics
    
    The function is optimized for:
    - Memory efficiency: Processes file line by line
    - Speed: Uses chromosome-based indexing for region lookup
    - Accuracy: Prefers CIGAR-based coverage when available
    
    Args:
        paf_file: Path to the PAF file
        regions: List of BED regions to analyze
        min_mapq: Minimum mapping quality threshold (default: 0)
        
    Returns:
        Dictionary with processing statistics:
        - alignment_count: Total valid alignments processed
        - filtered_count: Alignments filtered by low quality
        - no_cigar_count: Alignments without CIGAR strings
        - processed_chroms: Set of chromosomes seen in alignments
    """
    # Organize regions by chromosome for efficient lookup
    # This avoids checking all regions for every alignment
    regions_by_chrom = defaultdict(list)
    for region in regions:
        regions_by_chrom[region.chrom].append(region)
    
    # Initialize processing statistics
    stats = {
        'alignment_count': 0,        # Total alignments processed
        'filtered_count': 0,         # Alignments filtered by low quality
        'no_cigar_count': 0,        # Alignments without CIGAR strings
        'processed_chroms': set()    # Unique chromosomes seen
    }
    
    try:
        with open(paf_file, 'r') as f:
            # Process file line by line for memory efficiency
            for line_num, line in enumerate(f, 1):
                # Skip empty lines
                if not line.strip():
                    continue
                
                # Parse PAF record
                paf_record = parse_paf_line(line)
                if not paf_record:
                    print(f"Warning: Line {line_num} has invalid format, skipping", 
                          file=sys.stderr)
                    continue
                
                # Get target chromosome
                target_chrom = paf_record['target_name']
                
                # Skip if no regions on this chromosome
                # This is a major optimization for large files
                if target_chrom not in regions_by_chrom:
                    continue
                
                # Track processed chromosomes
                stats['processed_chroms'].add(target_chrom)
                
                # Apply mapping quality filter
                # Low quality alignments can introduce noise in coverage
                if paf_record['mapq'] < min_mapq:
                    stats['filtered_count'] += 1
                    continue
                
                # Get regions on this chromosome
                chrom_regions = regions_by_chrom[target_chrom]
                
                # Calculate coverage based on CIGAR availability
                if 'cigar' not in paf_record:
                    # No CIGAR string - use simple block coverage
                    # This is less accurate as it assumes no gaps
                    stats['no_cigar_count'] += 1
                    
                    # Simple coverage: mark all bases in alignment block as covered
                    for region in chrom_regions:
                        # Calculate overlap between alignment and region
                        overlap_start = max(paf_record['target_start'], region.start)
                        overlap_end = min(paf_record['target_end'], region.end)
                        
                        if overlap_start < overlap_end:
                            # Update coverage for overlapping bases
                            for pos in range(overlap_start, overlap_end):
                                region.coverage[pos] += 1
                else:
                    # CIGAR-based coverage (more accurate)
                    # Parse CIGAR string and calculate precise coverage
                    cigar_ops = parse_cigar(paf_record['cigar'])
                    calculate_coverage_from_cigar(
                        paf_record['target_start'], 
                        cigar_ops, 
                        chrom_regions,
                        target_chrom
                    )
                
                # Update alignment count
                stats['alignment_count'] += 1
                
                # Progress reporting for large files
                if line_num % 100000 == 0:
                    print(f"Processed {line_num:,} lines, {stats['alignment_count']:,} alignments", 
                          file=sys.stderr)
                
    except FileNotFoundError:
        print(f"Error: File not found - {paf_file}", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error: Failed to read file - {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nProcessing interrupted by user", file=sys.stderr)
        sys.exit(1)
    
    return stats


def print_summary(title: str, stats_list: List[Dict], process_stats: Dict, 
                 bed1_count: int = None, bed2_count: int = None):
    """
    Print summary statistics for a set of regions.
    
    This function generates a comprehensive summary report that provides
    an overview of coverage across all analyzed regions. It includes both
    processing statistics and coverage metrics.
    
    The summary helps answer high-level questions:
    - How many regions were analyzed?
    - What percentage of regions are fully covered?
    - What is the overall coverage breadth and depth?
    - How many alignments contributed to the coverage?
    
    Args:
        title: Title for the summary section
        stats_list: List of region statistics from calculate_region_stats()
        process_stats: Processing statistics from process_paf_file()
        bed1_count: Number of regions from first BED file (optional)
        bed2_count: Number of regions from second BED file (optional)
    """
    print(f"\n=== {title} Coverage Summary ===")
    
    # Show BED file breakdown if analyzing multiple files
    if bed1_count is not None and bed2_count is not None:
        print(f"BED1 regions: {bed1_count}")
        print(f"BED2 regions: {bed2_count}")
    
    # Basic counts
    print(f"Total regions: {len(stats_list)}")
    print(f"Alignments processed: {process_stats['alignment_count']:,}")
    print(f"Low quality filtered: {process_stats['filtered_count']:,}")
    print(f"No CIGAR information: {process_stats['no_cigar_count']:,}")
    
    # Show which chromosomes were involved
    if process_stats['processed_chroms']:
        chroms = sorted(process_stats['processed_chroms'])
        if len(chroms) <= 10:
            print(f"Chromosomes involved: {', '.join(chroms)}")
        else:
            print(f"Chromosomes involved: {len(chroms)} chromosomes")
    
    # Calculate overall coverage statistics
    total_length = sum(s['region'].length for s in stats_list)
    total_covered = sum(s['covered_bases'] for s in stats_list)
    
    # Categorize regions by coverage level
    fully_covered = sum(1 for s in stats_list if s['coverage_percent'] == 100)
    partially_covered = sum(1 for s in stats_list if 0 < s['coverage_percent'] < 100)
    not_covered = sum(1 for s in stats_list if s['coverage_percent'] == 0)
    
    # Overall coverage metrics
    print(f"\nTotal length: {total_length:,} bp")
    print(f"Covered bases: {total_covered:,} bp")
    if total_length > 0:
        print(f"Overall coverage: {total_covered/total_length*100:.2f}%")
    
    # Region coverage breakdown
    print(f"\nRegion coverage breakdown:")
    print(f"  Fully covered (100%): {fully_covered} regions")
    print(f"  Partially covered (0-100%): {partially_covered} regions")
    print(f"  Not covered (0%): {not_covered} regions")
    
    # Depth statistics across all regions
    avg_depths = [s['average_depth'] for s in stats_list if s['average_depth'] > 0]
    if avg_depths:
        print(f"\nAverage depth statistics:")
        print(f"  Mean: {sum(avg_depths)/len(avg_depths):.2f}×")
        print(f"  Min: {min(avg_depths):.2f}×")
        print(f"  Max: {max(avg_depths):.2f}×")
        
        # Additional depth percentiles could be calculated here
        # For example: median, quartiles, etc.


def print_region_stats(stats: Dict, source: str = ""):
    """
    Print detailed statistics for a single region.
    
    This function provides a detailed view of coverage for an individual
    region, which is useful for:
    - Identifying specific problematic regions
    - Understanding coverage patterns
    - Quality control and validation
    
    Args:
        stats: Region statistics from calculate_region_stats()
        source: Optional source label (e.g., "BED1", "BED2")
    """
    region = stats['region']
    
    # Format source label if provided
    prefix = f"[{source}] " if source else ""
    
    # Region identification
    print(f"\n{prefix}{region.name}")
    
    # Region coordinates (displayed in 1-based for readability)
    print(f"  Location: {region}")
    
    # Region size
    print(f"  Length: {region.length:,} bp")
    
    # Coverage breadth
    print(f"  Coverage: {stats['covered_bases']:,} bp ({stats['coverage_percent']:.2f}%)")
    
    # Coverage depth - overall
    print(f"  Average depth: {stats['average_depth']:.2f}×")
    
    # Coverage depth - covered bases only
    # This metric is useful when coverage is partial
    if stats['covered_bases'] > 0:
        print(f"  Average depth (covered only): {stats['average_depth_covered_only']:.2f}×")
    
    # Maximum depth observed
    print(f"  Maximum depth: {stats['max_depth']}×")
    
    # Additional statistics could be added here:
    # - Depth distribution histogram
    # - Coverage gaps analysis
    # - Depth variance/standard deviation


def main():
    """
    Main function - parse command line arguments and coordinate coverage analysis.
    
    This function orchestrates the entire analysis workflow:
    1. Parse and validate command line arguments
    2. Load BED regions from input files
    3. Process PAF alignments and calculate coverage
    4. Generate and display coverage statistics
    5. Output optional detailed reports
    
    The tool supports flexible analysis modes:
    - Separate analysis of two BED files
    - Merged analysis treating both BED files as one set
    - Detailed per-region reporting
    - Identification of uncovered regions
    """
    # Create argument parser with detailed help information
    parser = argparse.ArgumentParser(
        description='Calculate coverage statistics from PAF alignments for BED-defined regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Description:
  This tool calculates base-level coverage statistics from PAF alignment files
  for genomic regions defined in BED files. It supports accurate coverage
  calculation using CIGAR strings and can analyze multiple BED files separately
  or merged together.

Examples:
  # Basic coverage analysis for two BED files
  %(prog)s alignment.paf regions1.bed regions2.bed
  
  # Merge BED files and analyze as one set
  %(prog)s alignment.paf regions1.bed regions2.bed --merge-beds
  
  # Filter low quality alignments and show per-region details
  %(prog)s alignment.paf regions1.bed regions2.bed --min-mapq 30 --output-per-region
  
  # Find all uncovered regions
  %(prog)s alignment.paf regions1.bed regions2.bed --output-uncovered

Input formats:
  PAF: Pairwise Alignment Format from minimap2, minigraph, or similar tools
  BED: Browser Extensible Data format with 0-based, half-open intervals

Output:
  - Summary statistics for overall coverage
  - Optional per-region detailed statistics
  - Optional list of uncovered regions

Notes:
  - BED coordinates are 0-based, half-open intervals
  - PAF coordinates are also 0-based
  - CIGAR strings provide more accurate coverage than alignment blocks
  - Mapping quality filtering helps reduce noise from poor alignments
        """
    )
    
    # Positional arguments
    parser.add_argument('paf_file', 
                       help='PAF alignment file')
    parser.add_argument('bed_file1', 
                       help='First BED file with regions to analyze')
    parser.add_argument('bed_file2', 
                       help='Second BED file with regions to analyze')
    
    # Optional arguments
    parser.add_argument('--min-mapq', type=int, default=0, 
                       help='Minimum mapping quality threshold (default: 0, range: 0-255)')
    parser.add_argument('--output-per-region', action='store_true',
                       help='Output detailed statistics for each region')
    parser.add_argument('--output-uncovered', action='store_true',
                       help='Output list of regions with no coverage')
    parser.add_argument('--merge-beds', action='store_true',
                       help='Merge both BED files and analyze as a single set')
    
    # Parse command line arguments
    args = parser.parse_args()
    
    # Validate mapping quality parameter
    if args.min_mapq < 0 or args.min_mapq > 255:
        print(f"Error: --min-mapq must be between 0 and 255", file=sys.stderr)
        sys.exit(1)
    
    # Parse BED files
    print("Reading BED files...", file=sys.stderr)
    regions1 = parse_bed_file(args.bed_file1)
    regions2 = parse_bed_file(args.bed_file2)
    
    print(f"BED file 1: {len(regions1)} regions", file=sys.stderr)
    print(f"BED file 2: {len(regions2)} regions", file=sys.stderr)
    
    # Process based on merge option
    if args.merge_beds:
        # === Merged analysis mode ===
        # Treat both BED files as a single set of regions
        print("\nProcessing PAF file (merge mode)...", file=sys.stderr)
        
        # Combine all regions
        all_regions = regions1 + regions2
        
        # Process PAF file once for all regions
        process_stats = process_paf_file(args.paf_file, all_regions, args.min_mapq)
        
        # Calculate statistics for all regions
        print("\nCalculating coverage statistics...", file=sys.stderr)
        all_stats = [calculate_region_stats(region) for region in all_regions]
        
        # Output merged results
        print_summary("Merged Results", all_stats, process_stats, 
                     len(regions1), len(regions2))
        
        # Output per-region details if requested
        if args.output_per_region:
            print("\n=== Per-Region Statistics ===")
            # Track which BED file each region came from
            for i, stats in enumerate(all_stats):
                source = "BED1" if i < len(regions1) else "BED2"
                print_region_stats(stats, source)
                
    else:
        # === Separate analysis mode ===
        # Analyze each BED file independently
        
        # Process first BED file
        print("\nProcessing PAF file (BED1)...", file=sys.stderr)
        process_stats1 = process_paf_file(args.paf_file, regions1, args.min_mapq)
        stats1 = [calculate_region_stats(region) for region in regions1]
        
        # Process second BED file
        print("\nProcessing PAF file (BED2)...", file=sys.stderr)
        process_stats2 = process_paf_file(args.paf_file, regions2, args.min_mapq)
        stats2 = [calculate_region_stats(region) for region in regions2]
        
        # Output separate results
        print_summary("BED File 1", stats1, process_stats1)
        print_summary("BED File 2", stats2, process_stats2)
        
        # Output per-region details if requested
        if args.output_per_region:
            print("\n=== BED1 Per-Region Statistics ===")
            for stats in stats1:
                print_region_stats(stats, "BED1")
            
            print("\n=== BED2 Per-Region Statistics ===")
            for stats in stats2:
                print_region_stats(stats, "BED2")
    
    # Output uncovered regions if requested
    if args.output_uncovered:
        print("\n=== Uncovered Regions ===")
        
        # Collect statistics based on analysis mode
        if args.merge_beds:
            # Use merged statistics
            all_stats_combined = all_stats
        else:
            # Combine statistics from both BED files
            all_stats_combined = stats1 + stats2
        
        # Find regions with zero coverage
        uncovered_regions = []
        for stats in all_stats_combined:
            if stats['coverage_percent'] == 0:
                uncovered_regions.append(stats['region'])
        
        # Report uncovered regions
        if uncovered_regions:
            print(f"Found {len(uncovered_regions)} completely uncovered regions:")
            
            # Group by chromosome for better organization
            uncovered_by_chrom = defaultdict(list)
            for region in uncovered_regions:
                uncovered_by_chrom[region.chrom].append(region)
            
            # Display by chromosome
            for chrom in sorted(uncovered_by_chrom.keys()):
                print(f"\n  Chromosome {chrom}:")
                for region in uncovered_by_chrom[chrom]:
                    print(f"    {region.name}: {region}")
        else:
            print("All regions have at least partial coverage")
    
    print("\nAnalysis complete.", file=sys.stderr)


if __name__ == '__main__':
    # Script entry point
    # This ensures main() is only called when script is run directly,
    # not when imported as a module
    main()
