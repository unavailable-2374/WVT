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
   ``` 
   git clone https://github.com/unavailable-2374/WVT.git
   cd genome-analysis   
   ```
3. **Create a virtual environment** (recommended)  
   ```
   python3 -m venv venv  
   source venv/bin/activate  
   ```
4. **Install dependencies**  
   ```
   pip install -r requirements.txt  
   ```
---

## Usage

### PAF Dual Coverage

    python paf_dual_coverage.py \
      paf_file, query_bed, target_bed \
      -p threads

### Gene Cross-Matching Analysis

    python gene_cross_matching_analysis.py paf_file query_bed target_bed

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
