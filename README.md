# PAF Analysis Tools

A comprehensive toolkit for analyzing PAF (Pairwise Alignment Format) files with BED region coverage assessment and gene cross-matching analysis.

## Overview

This repository contains two powerful Python scripts for analyzing genomic alignments:

1. **`paf_dual_coverage.py`** - Evaluates PAF alignment coverage on both query and target BED regions
2. **`gene_cross_matching_analysis.py`** - Performs precise gene-to-gene matching analysis using CIGAR information

## Features

### PAF Dual Coverage Analysis (`paf_dual_coverage.py`)
- **Dual Coverage Assessment**: Simultaneously evaluates coverage on both query (BED1) and target (BED2) sequences
- **CIGAR-based Analysis**: Uses CIGAR strings for precise alignment coverage calculation
- **Multi-processing Support**: Parallel processing for large PAF files
- **Comprehensive Statistics**: Detailed coverage metrics including depth distribution and uncovered regions
- **Strand-aware**: Handles both forward and reverse strand alignments correctly
- **Quality Filtering**: Supports minimum mapping quality thresholds

### Gene Cross-matching Analysis (`gene_cross_matching_analysis.py`)
- **Precise Gene Mapping**: Maps query genes to target regions using CIGAR information
- **Cross-gene Analysis**: Identifies matches between different genes
- **Quality Thresholds**: Configurable gene coverage and identity thresholds
- **Diversity Analysis**: Analyzes gene matching patterns and diversity
- **Detailed Output**: Comprehensive matching statistics and alignment details
- **Debug Support**: Built-in debugging for specific genes

## Installation

### Requirements
- Python 3.6 or higher
- Required packages: `tqdm`, `argparse`, `collections`, `concurrent.futures`

### Setup
```bash
git clone https://github.com/unavailable-2374/WVT
cd WVT
chmod +x *.py
```

## Usage

### PAF Dual Coverage Analysis

#### Basic Usage
```bash
./paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed
```

#### Advanced Usage
```bash
# Use 40 processes for parallel processing
./paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed -p 40

# Filter low-quality alignments and output detailed statistics
./paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed \
    --min-mapq 30 --output-per-region --output-uncovered

# Save results to file
./paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed \
    -o coverage_results.txt --output-per-region
```

#### Parameters
- `paf_file`: PAF alignment file path
- `query_bed`: BED file for query sequences (BED1)
- `target_bed`: BED file for target sequences (BED2)
- `--min-mapq`: Minimum mapping quality threshold (default: 0)
- `--processes, -p`: Number of processes (default: CPU count)
- `--output-per-region`: Output detailed statistics for each region
- `--output-uncovered`: Output detailed information about uncovered regions
- `--output-file, -o`: Output file path (default: stdout)

### Gene Cross-matching Analysis

#### Basic Usage
```bash
./gene_cross_matching_analysis.py alignment.paf query_genes.bed target_genes.bed
```

#### Advanced Usage
```bash
# Set custom thresholds
./gene_cross_matching_analysis.py alignment.paf query_genes.bed target_genes.bed \
    --min-gene-coverage 0.8 --min-gene-identity 0.9

# Generate output files
./gene_cross_matching_analysis.py alignment.paf query_genes.bed target_genes.bed \
    --output detailed_matches.txt --matrix gene_matrix.txt

# Debug specific gene
./gene_cross_matching_analysis.py alignment.paf query_genes.bed target_genes.bed \
    --debug-gene GENE_NAME --show-gene-details GENE_NAME
```

#### Parameters
- `paf_file`: PAF alignment file with CIGAR information
- `query_bed`: BED file containing query gene regions
- `target_bed`: BED file containing target gene regions
- `--min-gene-coverage`: Minimum gene coverage threshold (default: 0.8)
- `--min-gene-identity`: Minimum gene identity threshold (default: 0.9)
- `--output`: Output detailed matching information to file
- `--matrix`: Output gene matching matrix to file
- `--debug-gene`: Debug specific gene matching
- `--show-gene-details`: Show detailed matching information for specific gene
- `--show-alignment`: Show detailed mapping for specific PAF line (provide line number)

## Input File Formats

### PAF Files
PAF files should contain standard PAF format with CIGAR strings in the optional fields:
```
query_name    query_length    query_start    query_end    strand    target_name    target_length    target_start    target_end    matches    alignment_length    mapping_quality    cg:Z:CIGAR_STRING
```

### BED Files
Standard BED format with at least 4 columns:
```
chromosome    start    end    name
```

For gene analysis, gene names are extracted from the name field using the pattern `prefix_genename`.

## Output Formats

### Coverage Analysis Output
```
=== Query Sequence Coverage Statistics (BED1) ===
Region count: 100
Processed alignments: 50,000
Low quality filtered: 1,000
No CIGAR information: 500

Total length: 1,000,000 bp
Covered bases: 950,000 bp
Overall coverage: 95.00%

Region coverage:
  Fully covered (100%): 80 regions
  Partially covered (0-100%): 15 regions
  Uncovered (0%): 5 regions
```

### Gene Matching Output
```
=== Gene Matching Type Statistics ===
Same gene matches: 150 pairs
Cross-gene matches: 25 pairs
Total matches: 175 pairs

Cross-gene Match Details (sorted by best identity):
Query Gene               Target Gene              Matches   Total Length   Avg Identity   Best Identity
GENE_A                   GENE_B                        3         1,500         92.5%         95.2%
GENE_C                   GENE_D                        2         1,200         89.8%         91.5%
```

## Algorithm Details

### Coverage Calculation
1. **Interval Tree Construction**: BED regions are indexed using interval trees for efficient overlap detection
2. **CIGAR Processing**: Each alignment is processed using CIGAR operations to determine exact coverage
3. **Strand Handling**: Reverse strand alignments are correctly mapped to forward coordinates
4. **Parallel Processing**: Large PAF files are split into chunks for parallel processing

### Gene Matching Algorithm
1. **CIGAR-based Mapping**: Query gene positions are mapped to target coordinates using CIGAR operations
2. **Overlap Detection**: Target genes overlapping with mapped regions are identified
3. **Quality Assessment**: Gene coverage and identity are calculated within gene boundaries
4. **Threshold Filtering**: Only matches meeting coverage and identity thresholds are reported

## Performance Considerations

### Memory Usage
- Memory usage scales with the number of BED regions and alignment density
- Large files benefit from increased process count for parallel processing
- Interval trees provide O(log n) lookup performance

### Processing Speed
- Typical processing speed: 10,000-50,000 alignments/second depending on CIGAR complexity
- Multi-processing provides near-linear speedup for large files
- CIGAR parsing is the main computational bottleneck

## Troubleshooting

### Common Issues

1. **Missing CIGAR Information**
   - Some alignments may lack CIGAR strings
   - The tool falls back to simple start-end coverage calculation
   - Use alignment tools that generate CIGAR strings (e.g., minimap2 with `-c` flag)

2. **Memory Issues**
   - Reduce the number of processes if memory is limited
   - Consider splitting large PAF files into smaller chunks

3. **Performance Issues**
   - Increase process count for CPU-bound operations
   - Ensure sufficient RAM for parallel processing

### Debug Options
```bash
# Debug specific gene matching
./gene_cross_matching_analysis.py alignment.paf query.bed target.bed --debug-gene GENE_NAME

# Show detailed alignment mapping
./gene_cross_matching_analysis.py alignment.paf query.bed target.bed --show-alignment 12345
```

## Example Workflows

### Workflow 1: Basic Coverage Analysis
```bash
# Step 1: Run coverage analysis
./paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed \
    --min-mapq 20 --output-per-region -o coverage_results.txt

# Step 2: Identify uncovered regions
./paf_dual_coverage.py alignment.paf query_regions.bed target_regions.bed \
    --output-uncovered -o uncovered_regions.txt
```

### Workflow 2: Gene Matching Analysis
```bash
# Step 1: Basic gene matching
./gene_cross_matching_analysis.py alignment.paf query_genes.bed target_genes.bed \
    --min-gene-coverage 0.8 --min-gene-identity 0.9

# Step 2: Generate detailed output
./gene_cross_matching_analysis.py alignment.paf query_genes.bed target_genes.bed \
    --output detailed_matches.txt --matrix gene_matrix.txt

# Step 3: Debug specific problematic gene
./gene_cross_matching_analysis.py alignment.paf query_genes.bed target_genes.bed \
    --debug-gene PROBLEM_GENE --show-gene-details PROBLEM_GENE
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use these tools in your research, please cite:
```
PAF Analysis Tools: Comprehensive toolkit for PAF alignment analysis with BED region coverage assessment
```

## Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Check the troubleshooting section
- Review the example workflows

## Changelog

### Version 1.0.0
- Initial release with dual coverage analysis
- Gene cross-matching with CIGAR support
- Multi-processing capabilities
- Comprehensive output formats
