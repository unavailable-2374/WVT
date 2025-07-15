# PAF Dual Coverage & Gene Cross-Matching Analysis

This repository contains two complementary Python scripts for genomic alignment analysis:

1. **`paf_dual_coverage.py`**  
   - Reads PAF-format alignment files  
   - Computes “dual” coverage (both query- and reference-based coverage)  
   - Merges overlapping alignment intervals  
   - Outputs per-sequence coverage summaries  

2. **`gene_cross_matching_analysis.py`**  
   - Aggregates cross-species gene match data from TSV files  
   - Filters alignments by score and coverage thresholds  
   - Computes per-gene summary statistics (mean, median, max scores, match counts)  
   - Produces CSV summary and publication-ready plots  

---

## Table of Contents

- [Features](#features)  
- [Installation](#installation)  
- [Usage](#usage)  
  - [PAF Dual Coverage](#paf-dual-coverage)  
  - [Gene Cross-Matching Analysis](#gene-cross-matching-analysis)  
- [Command-Line Arguments](#command-line-arguments)  
- [Output Files](#output-files)  
- [Dependencies](#dependencies)  
- [Examples](#examples)  
- [Contributing](#contributing)  
- [License](#license)  

---

## Features

- **High-performance PAF parsing** using pure Python  
- **Efficient interval merging**, avoids double-counting overlaps  
- **Flexible filtering** by score & coverage in cross-species matching  
- **Data aggregation** via pandas for per-gene statistics  
- **Automated plotting** of distribution histograms and scatter plots  
- Well-documented code with extensive inline comments for easy customization  

---

## Installation

1. **Clone the repo**  
    
    git clone [https://github.com/your-org/genome-analysis.git ](https://github.com/unavailable-2374/WVT.git) 
    cd genome-analysis  

2. **Create a virtual environment** (recommended)  
    
    python3 -m venv venv  
    source venv/bin/activate  

3. **Install dependencies**  
    
    pip install -r requirements.txt  

---

## Usage

### PAF Dual Coverage

    python paf_dual_coverage.py \
      --paf path/to/alignments.paf \
      --out dual_coverage_summary.tsv

- `--paf` — Input PAF file (tab-delimited, standard 12-column format)  
- `--out` — Output TSV; columns: `Name`, `Type` (`query` / `reference`), `Covered_bases`

### Gene Cross-Matching Analysis

From within a directory containing one or more `*.tsv` cross-match files:

    python gene_cross_matching_analysis.py

This will:

1. Read all `*.tsv` files in the current directory  
2. Filter matches to `score >= 50` and `coverage >= 0.8`  
3. Generate `gene_match_summary.csv` with per-gene metrics  
4. Save two PNG plots:  
   - `gene_match_summary_mean_score_histogram.png`  
   - `gene_match_summary_count_vs_score.png`

---

## Command-Line Arguments

### `paf_dual_coverage.py`

| Argument | Description                    | Required |
| -------- | ------------------------------ | -------- |
| `--paf`  | Path to input PAF file         | Yes      |
| `--out`  | Path to write coverage summary | Yes      |

### `gene_cross_matching_analysis.py`

No arguments. The script auto-detects `*.tsv` in the working directory.

---

## Output Files

- **`dual_coverage_summary.tsv`**  
      
      #Name    Type       Covered_bases  
      seq1     query      123456  
      seq1     reference  98765  
      ...  

- **`gene_match_summary.csv`**  
      
      query_gene,mean_score,median_score,max_score,match_count  
      GeneA,75.3,74.0,98.1,5  
      GeneB,62.1,60.5,85.0,3  
      ...  

- **Plots**  
  - `*_mean_score_histogram.png` — histogram of per-gene mean scores  
  - `*_count_vs_score.png`    — scatter of match count vs. mean score  

---

## Dependencies

- Python ≥ 3.6  
- `numpy`  
- `pandas`  
- `matplotlib`  
- `argparse` (stdlib)  
- `glob` (stdlib)  

Optionally, via conda:

    conda create -n genome-analysis python=3.9 numpy pandas matplotlib  
    conda activate genome-analysis  

---

## Examples

1. **Compute dual coverage** for a set of long-read alignments:  
      
      python paf_dual_coverage.py \
        --paf human_chr1_alignments.paf \
        --out chr1_dual_coverage.tsv  

2. **Analyze gene matches** across rice and grape species:  
      
      mv rice_vs_grape.tsv .  
      mv grape_vs_maize.tsv .  
      python gene_cross_matching_analysis.py  

3. **Customize thresholds**  
   Edit the `min_score` and `min_coverage` defaults in  
   `gene_cross_matching_analysis.py` to suit your data quality.

---

## Contributing

1. Fork this repository  
2. Create a feature branch (`git checkout -b feature/xyz`)  
3. Commit your changes (`git commit -m "Add feature"`)  
4. Push to your fork (`git push origin feature/xyz`)  
5. Open a Pull Request against `main`  

Please adhere to the existing coding style and include tests for any new functionality.

---

## License

Released under the [MIT License](LICENSE).  
Feel free to use, modify, and distribute as you see fit.
