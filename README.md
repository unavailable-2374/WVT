# PAF Coverage Tools

A collection of Python tools for analyzing alignment coverage from PAF (Pairwise Alignment Format) files against BED-defined regions. These tools provide comprehensive base-level coverage statistics, supporting both query and target sequence analysis with multi-core processing capabilities.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Tools](#tools)
  - [1. Query-Target Coverage Calculator](#1-query-target-coverage-calculator-1py)
  - [2. General Coverage Calculator](#2-general-coverage-calculator-coveragepy)
- [Input Formats](#input-formats)
- [Usage Examples](#usage-examples)
- [Output Format](#output-format)
- [Performance Considerations](#performance-considerations)
- [Contributing](#contributing)
- [License](#license)

## Overview

This repository contains two specialized tools for analyzing sequence alignment coverage:

1. **`paf_dual_coverage.py`** - A high-performance tool specifically designed for analyzing query and target coverage from PAF alignments with multi-core processing support
2. **`paf_coverage.py`** - A general-purpose coverage calculator for any BED-defined regions

Both tools parse CIGAR strings for accurate base-level coverage calculation and provide detailed statistics including coverage percentage, depth distribution, and per-region analysis.

## Features

### Common Features
- Base-level coverage calculation using CIGAR string parsing
- Support for PAF format (minimap2, minigraph, etc.)
- BED file support for defining regions of interest
- Mapping quality filtering
- Detailed coverage statistics (coverage %, average depth, max depth)
- Per-region detailed reports
- Identification of uncovered regions

### Tool-Specific Features

**Query-Target Coverage Calculator (`paf_dual_coverage.py`)**:
- Multi-core parallel processing for large PAF files
- Separate tracking of query and target coverage
- Interval tree implementation for efficient region lookup
- Progress bar with real-time statistics
- Optimized for high-throughput analysis

**General Coverage Calculator (`paf_coverage.py`)**:
- Support for merging multiple BED files
- Flexible coverage analysis for any reference sequences
- Simpler single-threaded implementation

## Requirements

- Python 3.6 or higher
- Required Python packages:
  - `tqdm` (for progress bars in `paf_dual_coverage.py`)
  - Standard library modules: `argparse`, `re`, `sys`, `os`, `collections`, `typing`, `concurrent.futures`, `multiprocessing`, `time`, `bisect`, `dataclasses`

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/paf-coverage-tools.git
cd paf-coverage-tools
```

2. Install dependencies:
```bash
pip install tqdm
```

3. Make scripts executable:
```bash
chmod +x paf_dual_coverage.py paf_coverage.py
```

## Tools

### 1. Query-Target Coverage Calculator (`paf_dual_coverage.py`)

This tool is designed for analyzing alignments where you need to track coverage on both query and target sequences separately.

**Basic Usage:**
```bash
./paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed
```

**Options:**
- `--min-mapq`: Minimum mapping quality (default: 0)
- `--processes, -p`: Number of processes for parallel processing (default: CPU count)
- `--output-per-region`: Output detailed statistics for each region
- `--output-uncovered`: Output detailed information about uncovered regions
- `--output-file, -o`: Output file path (default: stdout)

### 2. General Coverage Calculator (`paf_coverage.py`)

This tool provides flexible coverage analysis for any BED-defined regions against a reference.

**Basic Usage:**
```bash
./paf_coverage.py alignment.paf regions1.bed regions2.bed
```

**Options:**
- `--min-mapq`: Minimum mapping quality (default: 0)
- `--output-per-region`: Output detailed statistics for each region
- `--output-uncovered`: Output detailed information about uncovered regions
- `--merge-beds`: Merge statistics from both BED files

## Input Formats

### PAF Format
Standard PAF format as produced by minimap2, wfmash, or other aligners:
```
query_name query_len query_start query_end strand target_name target_len target_start target_end matches alignment_len mapq [optional_fields]
```

The tools specifically look for CIGAR strings in the optional fields (format: `cg:Z:cigar_string`).

### BED Format
Standard BED format (0-based, half-open intervals):
```
chromosome start end [name]
```

Example:
```
chr1    1000    2000    region1
chr1    5000    6000    region2
chr2    1000    3000    region3
```

## Usage Examples

### Example 1: Basic Query-Target Analysis
```bash
# Analyze coverage with default settings
./paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed

# Use 40 cores for processing large files
./paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed -p 40

# Filter low-quality alignments and save detailed output
./paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed \
    --min-mapq 30 \
    --output-per-region \
    --output-file results.txt
```

### Example 2: General Coverage Analysis
```bash
# Basic coverage analysis
./paf_coverage.py alignment.paf regions1.bed regions2.bed

# Merge BED files and analyze together
./paf_coverage.py alignment.paf regions1.bed regions2.bed --merge-beds

# Output all uncovered regions
./paf_coverage.py alignment.paf regions1.bed regions2.bed \
    --output-uncovered \
    --min-mapq 20
```

### Example 3: Pipeline Integration
```bash
# Run minimap2 and analyze coverage in a pipeline
wfmash -t 8 reference.fa query.fa | \
    ./paf_dual_coverage.py - query_regions.bed target_regions.bed -p 8
```

## Output Format

### Summary Statistics
Both tools provide comprehensive summary statistics:
- Total number of regions
- Number of processed alignments
- Filtered alignments (low quality)
- Alignments without CIGAR information
- Overall coverage percentage
- Coverage distribution (fully covered, partially covered, uncovered regions)
- Average depth statistics

### Per-Region Statistics (with `--output-per-region`)
- Region name and coordinates
- Region length
- Covered bases and percentage
- Average depth (overall and covered bases only)
- Maximum depth

### Example Output
```
=== Query Coverage Statistics (BED1) ===
Regions: 150
Processed alignments: 10,543
Low quality filtered: 523
No CIGAR info: 0
Query sequences involved: 5

Total length: 1,500,000 bp
Covered bases: 1,350,000 bp
Overall coverage: 90.00%

Region coverage:
  Fully covered (100%): 120 regions
  Partially covered (0-100%): 25 regions
  Uncovered (0%): 5 regions

Average depth statistics:
  Average: 3.45×
  Minimum: 0.50×
  Maximum: 12.30×
```

## Performance Considerations

### Memory Usage
- Both tools use interval trees for efficient region lookup
- Memory usage scales with the number of regions in BED files
- The parallel version (`paf_dual_coverage.py`) duplicates data structures across processes

### Processing Speed
- `paf_dual_coverage.py` can process millions of alignments per minute on modern multi-core systems
- Performance scales nearly linearly with the number of CPU cores
- I/O typically becomes the bottleneck for very large files

### Optimization Tips
1. Use more processes (`-p`) for large PAF files
2. Pre-filter PAF files by mapping quality if not all alignments are needed
3. Consider splitting very large BED files if memory is limited
4. Use SSDs for better I/O performance with large files

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

### Development Guidelines
1. Follow PEP 8 style guidelines
2. Add docstrings to all functions
3. Include type hints where appropriate
4. Update documentation for new features
5. Add tests for new functionality

## Contact

For questions, issues, or suggestions, please open an issue on GitHub or contact [scao7@uthsc.edu].
