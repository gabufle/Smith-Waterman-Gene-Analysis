# Smith-Waterman Gene Analysis

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Topic: Bioinformatics](https://img.shields.io/badge/Topic-Bioinformatics-green.svg)](https://github.com/topics/bioinformatics)

An implementation and analytical toolkit for local sequence alignment of biological sequences (DNA, RNA, and proteins) using the **Smith-Waterman Algorithm**. This repository provides tools for dynamic programming matrix computation, traceback path reconstruction, scoring customization, and comparative gene analysis.

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [How It Works](#how-it-works)
- [Installation](#installation)
- [Usage](#usage)
  - [Command Line Interface](#command-line-interface)
  - [Python API](#python-api)
- [Configuration & Scoring](#configuration--scoring)
- [Example Output](#example-output)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)

---

## Overview

The **Smith-Waterman algorithm** is a classic dynamic programming approach used in bioinformatics to find the optimal local alignment between two nucleotide or amino acid sequences. Unlike global alignment algorithms (e.g., Needleman-Wunsch), Smith-Waterman identifies highly conserved local regions between sequences of varying lengths or high divergence.

This repository provides a flexible Python framework to perform local sequence alignments on FASTA files or raw gene strings, calculate optimal alignment scores, visualize tracebacks, and evaluate sequence similarity.

---

## Key Features

- **Local Alignment Core**: Full dynamic programming implementation of the Smith-Waterman algorithm.
- **Flexible Input Formats**: Load sequences from `.fasta`, `.txt`, or directly via Python objects.
- **Customizable Scoring Rules**:
  - Support for custom Match / Mismatch scores.
  - Linear and Affine gap penalty models (Gap Open and Gap Extension).
  - Built-in substitution matrices for protein sequences (e.g., `BLOSUM62`, `PAM250`).
- **Traceback Reconstruction**: Automatically tracks and reconstructs single or multiple optimal local alignment paths.
- **Alignment Metrics**: Outputs total alignment score, sequence identity percentage, sequence similarity, and alignment visualization.
- **Matrix Visualization**: (Optional) Render dynamic programming scoring heatmaps and traceback paths.

---

## How It Works

The algorithm constructs a dynamic programming matrix $H$ of size $(m+1) \times (n+1)$, where $m$ and $n$ are the lengths of sequence $A$ and sequence $B$, respectively.

### Recurrence Relation

For each cell $H(i, j)$, the score is calculated as:

$$H(i, j) = \max \begin{cases} 0 \\ H(i-1, j-1) + s(A_i, B_j) & \text{(Match / Mismatch)} \\ \max_{k \ge 1} \{ H(i-k, j) - W_k \} & \text{(Deletion in } B \text{)} \\ \max_{l \ge 1} \{ H(i, j-l) - W_l \} & \text{(Insertion in } B \text{)} \end{cases}$$

Where:
- $s(A_i, B_j)$ represents the substitution score for aligning characters $A_i$ and $B_j$.
- $W_k$ is the gap penalty function for a gap of length $k$.
- A minimum value of $0$ ensures that only locally optimal regions with positive scores are considered.

### Traceback

Traceback begins at the cell with the maximum value $\max(H(i, j))$ in the matrix and continues until reaching a cell with a score of $0$.

---

## Installation

### Prerequisites

- **Python**: Version 3.8 or higher.
- Required Python libraries (if visualizer/matrices are enabled): `numpy`, `matplotlib`, `biopython` (optional).

### Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/gabufle/Smith-Waterman-Gene-Analysis.git
   cd Smith-Waterman-Gene-Analysis
   ```

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

---

## Usage

### Command Line Interface

Run local alignment between two FASTA files:

```bash
python main.py --seq1 data/gene_a.fasta --seq2 data/gene_b.fasta --match 5 --mismatch -4 --gap 6
```

#### Command-Line Arguments

| Parameter | Short | Type | Default | Description |
|---|---|---|---|---|
| `--seq1` | `-1` | File / String | Required | Path to sequence 1 (FASTA/text) or raw sequence |
| `--seq2` | `-2` | File / String | Required | Path to sequence 2 (FASTA/text) or raw sequence |
| `--match` | `-m` | Int | `2` | Score for matching characters |
| `--mismatch` | `-x` | Int | `-1` | Penalty for mismatched characters |
| `--gap` | `-g` | Int | `-2` | Linear gap penalty |
| `--matrix` | `-mat` | String | `None` | Substitution matrix name (e.g., `BLOSUM62`) |
| `--output` | `-o` | File | `None` | Save alignment results to specified text file |

---

### Python API

You can also import and use the aligner within your own scripts:

```python
from src.smith_waterman import SmithWaterman

# Define sequences
seq_a = "HEAGAWGHEE"
seq_b = "PAWHEAE"

# Initialize aligner with linear gap penalty
aligner = SmithWaterman(
    match_score=3,
    mismatch_penalty=-1,
    gap_penalty=-2,
    matrix="BLOSUM62" # Optional for amino acid alignments
)

# Run alignment
alignment_result = aligner.align(seq_a, seq_b)

# Display results
print(f"Optimal Score: {alignment_result.score}")
print("Aligned Sequence A:", alignment_result.aligned_seq1)
print("                    ", alignment_result.match_line)
print("Aligned Sequence B:", alignment_result.aligned_seq2)
print(f"Identity: {alignment_result.identity_percent:.2f}%")
```

---

## Configuration & Scoring

Default scoring parameters can be modified in `config.py` or passed dynamically via function calls.

| Sequence Type | Parameter | Recommended Default |
|---|---|---|
| **DNA / RNA** | Match / Mismatch | $+5 / -4$ or $+2 / -3$ |
| **DNA / RNA** | Gap Penalty | $-5$ (Linear) or Open $-10$, Extend $-1$ (Affine) |
| **Protein** | Substitution Matrix | `BLOSUM62` or `PAM250` |
| **Protein** | Gap Penalty | Open $-11$, Extend $-1$ |

---

## Example Output

Executing an alignment on two sample gene fragments yields output formatted as follows:

```text
==================================================
        SMITH-WATERMAN LOCAL ALIGNMENT
==================================================
Sequence 1 Length: 28 bp
Sequence 2 Length: 25 bp
Scoring: Match = +5, Mismatch = -4, Gap = -6

Optimal Alignment Score: 45
Alignment Length: 10
Identity: 9 / 10 (90.0%)

Alignment Visualization:
Seq1:  12  G A C G G A T T A G  21
           | | | | | | |   |
Seq2:   8  G A C G G A T - A G  16
==================================================
```

---

## Project Structure

```text
Smith-Waterman-Gene-Analysis/
├── data/                  # Sample FASTA files and gene sequence data
│   ├── gene_sample_1.fasta
│   └── gene_sample_2.fasta
├── src/                   # Core source code
│   ├── __init__.py
│   ├── smith_waterman.py  # SW algorithm matrix generation and traceback
│   ├── fasta_parser.py    # Utilities for loading FASTA format sequences
│   └── utils.py           # Scoring matrices and metric calculations
├── tests/                 # Unit tests
│   └── test_alignment.py
├── config.py              # Configuration defaults
├── main.py                # CLI entry point
├── requirements.txt       # Project dependencies
└── README.md              # Project documentation
```

---

## Contributing

Contributions are welcome! If you find a bug or would like to add new feature functionality (such as parallelized SW matrix evaluation or additional visualization tools):

1. Fork the repository.
2. Create your feature branch (`git checkout -b feature/AmazingFeature`).
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`).
4. Push to the branch (`git checkout push origin feature/AmazingFeature`).
5. Open a Pull Request.

---

## License

Distributed under the MIT License. See `LICENSE` for more information.
