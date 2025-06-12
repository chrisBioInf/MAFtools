# MAFtools

A comprehensive Python toolkit for processing and analyzing Multiple Alignment Format (MAF) files.

## Overview

MAFtools is a command-line utility that provides various operations for MAF alignment files, including filtering, merging, windowing, and format conversion. It's particularly useful for comparative genomics analyses where you need to manipulate whole-genome alignments.

**Version:** 1.0.0  
**Author:** Christopher Klapproth  
**Institution:** University Leipzig  
**License:** GPLv2

### Requirements
- Python 3.x
- BioPython
- NumPy

```bash
pip install biopython numpy
```

## Usage

```bash
python MAFtools.py [program] [options] [MAF file]
```

Where `program` is one of: `describe`, `merge`, `filter`, `window`, `select`, `toBed`

## Core Functionality

### 1. describe
Provides summary statistics about a MAF file.

```bash
python MAFtools.py describe alignment.maf
```

**Output includes:**
- Number of unique identifiers and sequences
- Total number of alignment blocks
- Total columns and average block length
- Average sequences per block
- Aligned reference nucleotides

### 2. select
Extracts alignment blocks based on BED coordinates.

```bash
python MAFtools.py select -b regions.bed -o extracted.maf alignment.maf
```

**Options:**
- `-b, --bed`: BED file with coordinates to extract (required)
- `-r, --right-flank`: Add nucleotides to the right side of each coordinate (default: 0)
- `-l, --left-flank`: Add nucleotides to the left side of each coordinate (default: 0)
- `-o, --output`: Output MAF file (default: output.maf)

**Example:** If a BED entry is `species1 1560 1960` and a MAF block starts at position 1500 with length 500, the program extracts the overlapping region (positions 1560-1960).

### 3. filter
Filters MAF alignment blocks based on various criteria.

```bash
python MAFtools.py filter -g 0.25 -n 3 -o filtered.maf alignment.maf
```

**Options:**
- `-g, --max-gaps`: Maximum gap fraction allowed (default: 0.25)
- `-a, --max-similarity`: Maximum pairwise similarity to reference (default: 1.0)
- `-m, --min-similarity`: Minimum pairwise similarity to reference (default: 0.5)
- `-l, --min-length`: Minimum block length in columns (default: 40)
- `-n, --min-seqs`: Minimum number of sequences required (default: 3)
- `-o, --output`: Output MAF file (default: output.maf)

**Filter criteria:**
- Average gap fraction across all sequences
- Pairwise sequence similarity to the reference
- Number of sequences in the block
- Alignment block length

### 4. window
Applies a sliding window approach to create overlapping smaller alignment blocks.

```bash
python MAFtools.py window -l 120 -s 40 -o windowed.maf alignment.maf
```

**Options:**
- `-l, --length`: Target length of cut alignment blocks (default: 120)
- `-s, --step`: Step size for sliding window (default: 40)
- `-o, --output`: Output MAF file (default: output.maf)

**Behavior:**
- Blocks shorter than the window length are kept as-is
- Longer blocks are cut into overlapping windows
- Example: A 200bp block with step=40 creates windows starting at positions 0, 40, 80, etc.

### 5. merge
Merges adjacent alignment blocks based on species consensus and distance criteria.

```bash
python MAFtools.py merge -s 0.75 -d 100 -o merged.maf alignment.maf
```

**Options:**
- `-s, --species-consensus`: Minimum consensus between blocks for merging (default: 0.75)
- `-r, --no-reference`: Flag if first sequence should NOT be considered reference
- `-d, --max-distance`: Maximum distance between blocks for merging (default: 0)
- `-l, --max-length`: Maximum length of merged blocks (default: 1000)
- `-g, --max-gaps`: Maximum gap fraction before dropping sequences (default: 0.5)
- `-m, --min-seqs`: Minimum shared sequences for merging (default: 2)
- `-o, --output`: Output MAF file (default: output.maf)

### 6. toBed
Converts MAF alignment blocks to BED format.

```bash
python MAFtools.py toBed -o alignment.bed alignment.maf
```

**Options:**
- `-o, --output_bed`: Output BED file (default: output.bed)

**Output format:**
Each line represents one alignment block's coordinates:
```
sequence_name    start    end
```

## Examples

### Extract specific genomic regions:

```
# Extract these regions from MAF
python MAFtools.py select -b Example/test.bed -o extracted.maf Example/test.maf
```

### Filter alignments for quality or similarity criteria:
```bash
# Keep blocks with <15% gaps, 3+ sequences, and 70-95% similarity
python MAFtools.py filter -g 0.15 -n 3 -m 0.70 -a 0.95 -o filtered.maf Example/test.maf
```

### Create sliding windows for analysis:
```bash
# Create 100bp windows with 25bp overlap
python MAFtools.py window -l 100 -s 25 -o windows.maf genome.maf
```

## MAF Format Notes

- Coordinates are 0-based, half-open intervals
- The reference sequence is typically the first sequence in each block
- Gap characters are represented by '-'
- Each block contains aligned sequences with metadata (start, size, strand, srcSize)

## Citation

If you use MAFtools in your research, please cite:
```
Klapproth, C. (2024). MAFtools: A toolkit for Multiple Alignment Format file manipulation. 
University Leipzig. Version 1.0.0.
```
