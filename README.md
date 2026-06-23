# ECSFinder

<p align="center">
    <img width="25%" src="https://user-images.githubusercontent.com/44384386/195381940-680064be-d53a-45b6-a5e1-a80ff1cb804e.jpg" alt="ECSFinder Logo"> 
</p>

<p align="center">
    <strong>A machine learning-powered tool for detecting evolutionarily conserved RNA secondary structures</strong>
</p>

---

## Overview

ECSFinder is a comprehensive bioinformatics pipeline that identifies and validates **evolutionarily conserved RNA secondary structures** (ECS) across multiple species alignments. By integrating structural prediction, statistical testing, and machine learning classification, ECSFinder provides robust detection and prioritization of functional RNA elements in genomic sequences.


### How It Works

```
MAF Files → Merge & Filter → Local Structure Prediction → Conservation Testing → 
Feature Extraction → Random Forest Classification → Validated ECS Elements
```


**Classification Features:**

The random forest model uses 9 key features:
- E-value from R-scape
- Number of significant covarying base pairs
- Minimum free energy (MFE)
- Pseudo-energy
- Mean pairwise identity (MPI)
- Standard deviation of null energies
- Average MFE from SISSIz background
- Z-score from SISSIz
- Structure conservation index (SCI)

---

## Table of Contents

- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Quick Start](#quick-start)
  - [Detailed Setup](#detailed-setup)
- [Usage](#usage)
  - [Basic Usage](#basic-usage)
  - [Command-Line Options](#command-line-options)
  - [Advanced Examples](#advanced-examples)
- [Reference Species Configuration](#reference-species-configuration)
- [Output Files](#output-files)
- [Interpreting Results](#interpreting-results)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Installation

### Prerequisites

**Required Software:**

| Tool | Version | Purpose |
|------|---------|---------|
| Java | 11+ | Run ECSFinder |
| SISSIz | 3.0 | Structural conservation testing |
| ViennaRNA | 2.5+ | RNA structure prediction |
| R-scape | Latest | Covariation analysis |
| R | 4.3+ | Random forest classification |
| Maven | 3.6+ | Build from source (optional) |

**R Packages:**
```r
install.packages(c("caret", "randomForest"))
```

### Quick Start

**Option 1: Use Prebuilt JAR (Recommended)**

```bash
# Clone the repository
git clone https://github.com/yourusername/ECSFinder.git
cd ECSFinder

# Run with the prebuilt JAR (using default Homo_sapiens reference)
java -jar target/ECSFinder-1.0.0.jar \
  -mafft
  -o output_directory \
  -i input.maf
```

**Option 2: Build from Source**

```bash
# Clone and build
git clone https://github.com/yourusername/ECSFinder.git
cd ECSFinder
mvn clean package


java -cp target/classes \
  ca.smithlab.vandalovejoy.ecsfinder.ECSFinder \
  [options] \
  -o output/directory \
  -i input.maf_or_dir

```

### Detailed Setup

#### 1. Install SISSIz 3.0
Use the directory 3_0_SISSIz and install SISSIz according to its README.

```bash

# Verify installation
which SISSIz
SISSIz --version
```


#### 2. Install ViennaRNA Package

```bash
# Download from https://www.tbi.univie.ac.at/RNA/
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.1.tar.gz

# Extract and install
tar -zxvf ViennaRNA-2.5.1.tar.gz
cd ViennaRNA-2.5.1
./configure
make
sudo make install

# Verify installation
which RNALalifold
which RNAalifold
RNALalifold --version
```

#### 3. Install R-scape

```bash
# Download from http://eddylab.org/R-scape/
# Install R-scape according to its README.

# Verify installation
which R-scape
R-scape -h
```

#### 4. Set Up R Environment

```bash
# Check R version (must be 4.3+)
R --version

# Install required packages
R -e 'install.packages(c("caret", "randomForest"), repos="http://cran.rstudio.com/")'

# Verify Rscript is available
which Rscript
```

#### 5. Verify ECSFinder Installation

```bash
cd ECSFinder

# Test with example data
java -jar target/ECSFinder-1.0.0.jar \
  -o test_output \
  -i example/test.maf \
  -c 2
```

---

## Usage

### Basic Usage

**Minimum required arguments:**

```bash
# Using default reference species (Homo_sapiens)
java -jar target/ECSFinder-1.0.0.jar \
  -o <output_directory> \
  -i <input.maf_or_directory>

# Or specify a different reference species
java -jar target/ECSFinder-1.0.0.jar \
  -ref <reference_species> \
  -o <output_directory> \
  -i <input.maf_or_directory>
```

**Example:**

```bash
# Using default Homo_sapiens reference
java -jar target/ECSFinder-1.0.0.jar \
  -o results/human_ECS \
  -i alignments/

# Using hg38 reference
java -jar target/ECSFinder-1.0.0.jar \
  -ref pan_paniscus \
  -o results/human_ECS \
  -i alignments/
```
Input requirements

Minimum sequences per alignment block: ECSFinder is designed for multi-species alignments and relies on evolutionary signal
for conserved RNA structure/covariation. Alignment blocks with <5 sequences are skipped or with sequences fewer than 50 nucleotides in length.

### Command-Line Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-i` | String | **Required** | Input MAF file or directory containing MAF files |
| `-o` | String | **Required** | Output directory for results |
| `-ref` | String | Homo_sapiens | Reference species identifier (see [Reference Species](#reference-species-configuration)) |
| `-c` | Integer | 4 | Number of CPU cores for parallel processing |
| `-g` | Integer | 50 | Maximum gap percentage allowed in alignments (0-100) |
| `-sszr` | Double | -3.0 | SISSIz Z-score threshold (more negative = more stringent) |
| `-mpi` | Integer | 50 | Minimum mean pairwise identity percentage (0-100) |
| `-mafft` | Flag | false | Realign each block with MAFFT-GINSI before analysis |
| `-v` | Flag | false | Enable verbose output for debugging |
| `-maxbp` | Integer | 300 | Maximum base pair span for RNALalifold |
| `-overlap` | Integer | 300 | Overlap length for block splitting |
| `-maxblock` | Integer | 5000 | Maximum block size for block splitting |
| `-threshold` | Double | 0.679 | RF model prediction threshold |

### Advanced Examples

#### Example 1: High-Stringency Analysis

Detect only the most conserved structures with high sequence identity:

```bash
java -jar target/ECSFinder-1.0.0.jar \
  -ref homo_sapiens \
  -sszr -5.0 \
  -mpi 70 \
  -g 30 \
  -c 16 \
  -o high_confidence_results \
  -i primate_alignments/
```

#### Example 2: Permissive Search

Cast a wider net for potential RNA structures:

```bash
java -jar target/ECSFinder-1.0.0.jar \
  -ref mm10 \
  -sszr -2.0 \
  -mpi 40 \
  -g 60 \
  -c 8 \
  -o exploratory_results \
  -i mammal_alignments/
```

#### Example 3: With Realignment

Use MAFFT to improve alignment quality before structure prediction:

```bash
java -jar target/ECSFinder-1.0.0.jar \
  -ref hg38 \
  -mafft \
  -c 12 \
  -o realigned_results \
  -i input.maf
```

#### Example 4: Verbose Debugging

Get detailed output for troubleshooting:

```bash
java -jar target/ECSFinder-1.0.0.jar \
  -ref hg38 \
  -v \
  -o debug_output \
  -i problem_alignment.maf \
  2>&1 | tee ecsfinder.log
```

---

## Reference Species Configuration

ECSFinder requires a **reference species** to convert alignment coordinates to genomic coordinates. By default, ECSFinder uses **Homo_sapiens** as the reference species. You can override this using the `-ref` parameter.

### Understanding MAF Species IDs

In MAF files, species are identified in `s` lines:

```
s Homo_sapiens.chr1  12345  200 + 248956422 ACUG...
s panTro6.chr1       23456  200 + 227471971 ACUG...
s rheMac10.chr1      34567  200 + 223616942 ACUG...
```

The species ID is everything before the first space (e.g., `Homo_sapiens.chr1`, `panTro6.chr1`).

**Default behavior:** If `-ref` is not specified, ECSFinder looks for sequences matching "Homo_sapiens".

### Substring Matching (Simple)

If `-ref` is a **plain string**, ECSFinder performs case-insensitive substring matching:

```bash
# Default - matches Homo_sapiens.chr1, Homo_sapiens.chr2, etc.
# (No -ref needed, this is the default)
java -jar target/ECSFinder-1.0.0.jar -o output -i input.maf

# Matches: hg38.chr1, hg38.chr2, hg38.chrX, etc.
-ref hg38

# Explicitly specify Homo_sapiens
-ref Homo_sapiens

# Matches: mm10.chr1, mm10.chr2, etc.
-ref mm10
```

### Regular Expression Matching (Advanced)

For precise control, wrap your pattern in `/slashes/` to use regex:

```bash
# Match only sequences starting with "homo_sapiens."
-ref "/^homo_sapiens\./"

# Match either hg38 or GRCh38 assemblies
-ref "/^(hg38|GRCh38)\./"

# Match human chromosome 1 only
-ref "/^hg38\.chr1$/"

# Match autosomes only (exclude X, Y, M)
-ref "/^hg38\.chr[0-9]+$/"
```

### Common Reference Species

| Species | Assembly | Substring | Regex Example |
|---------|----------|-----------|---------------|
| Human (default) | Homo_sapiens | `Homo_sapiens` | `/^Homo_sapiens\./` |
| Human | GRCh38/hg38 | `hg38` | `/^hg38\./` |
| Mouse | GRCm39/mm39 | `mm39` | `/^mm39\./` |
| Rat | mRatBN7 | `rn7` | `/^rn7\./` |
| Zebrafish | GRCz11 | `danRer11` | `/^danRer11\./` |
| Fruit fly | BDGP6 | `dm6` | `/^dm6\./` |

### Troubleshooting Reference Species

**Problem:** "No reference species found" warnings

```bash
# Use verbose mode to see what species IDs are in your MAF
java -jar target/ECSFinder-1.0.0.jar -v -o out -i input.maf 2>&1 | grep "species"

# Check your MAF file directly
grep "^s " input.maf | cut -d' ' -f2 | sort -u
```

**Solution:** Adjust `-ref` to match the actual species IDs in your file, or omit `-ref` to use the default Homo_sapiens reference.

---

## Output Files

ECSFinder creates several output files in the specified output directory:

### Directory Structure

```
output_directory/
├── output.maf              # Merged and filtered alignment
├── final.csv               # Feature matrix for all candidates
├── predictions.csv         # Random forest predictions
└── aln/                    # Alignment files (TP only)
    ├── chr1_12345_12400_10_85.2_0.2_0.45_55.1_0.3_450_+.aln
```

### Main Output Files

#### 1. `output.maf`

Merged MAF alignment with:
- Ancestor sequences removed
- Duplicate species filtered
- Ready for downstream analysis

#### 2. `final.csv`

Feature matrix with columns:

| Column | Description |
|--------|-------------|
| `name_file` | Alignment identifier |
| `min_energy` | Minimum free energy (kcal/mol) |
| `pseudo_energy` | Normalized pseudo-energy |
| `log_min_evalue` | Log₁₀ of minimum E-value from R-scape |
| `covarying_bp` | Number of significant covarying base pairs |
| `MPI` | Mean pairwise identity (%) |
| `average_MFE_sample` | Mean MFE from null background |
| `sd_sample` | Standard deviation of null energies |
| `zscore` | Z-score from SISSIz |
| `sci` | Structure conservation index |

#### 3. `predictions.csv`

Random forest predictions:

```csv
name_file,Predicted_Probabilities,Predicted_Class
chr5_RF00017_273_426_5_58.1_0.3_0.77_59.2_0.7_1317_+.aln,0.94,TP
chr5_RF00017_318_435_5_60.5_0.3_0.73_57.0_1.0_655_+.aln,0.84,TP
```

- **Predicted_Probabilities**: Confidence score (0-1)
- **Predicted_Class**: TP (true positive) or FP (false positive)

⚠️ **Note:** Only alignments classified as TP are retained in the `aln/` directory.

#### 4. Alignment Files (`aln/` directory)

For each true positive ECS:

**`.aln` file** - ClustalW format alignment

### Alignment Filename Convention

Filenames encode key features:

```
chr5_9958021_9958096_11_92.2_0.1_0.16_66.8_0.4_304_+.aln
│    │       │        │  │    │   │    │    │   │   │
│    │       │        │  │    │   │    │    │   │   └─ Strand
│    │       │        │  │    │   │    │    │   └───── Z-score × -100
│    │       │        │  │    │   │    │    └───────── Gap content
│    │       │        │  │    │   │    └────────────── GC content (%)
│    │       │        │  │    │   └─────────────────── Mean Shannon entropy
│    │       │        │  │    └─────────────────────── SD of MPI
│    │       │        │  └──────────────────────────── Mean pairwise identity (%)
│    │       │        └─────────────────────────────── Number of species
│    │       └──────────────────────────────────────── End coordinate
│    └──────────────────────────────────────────────── Start coordinate
└───────────────────────────────────────────────────── Chromosome
```

**Redirect to file:**

```bash
java -jar target/ECSFinder-1.0.0.jar ... > ecs_predictions.bed
```
---

## Interpreting Results


### Key Metrics


### Quality Control

**Good results typically show:**
- Z-scores < -3.0
- MPI between 50-90%
- Multiple covarying base pairs
- Probability > 0.7
- Low gap content (< 30%)

**Warning signs:**
- Very high MPI (> 95%) - may lack variation for covariation analysis
- Very low MPI (< 40%) - alignment may be unreliable
- High gap content (> 50%) - alignment quality issues
- No covarying pairs - weak structural evidence

---

## Troubleshooting

### Common Issues


**Problem: "SISSIz not found"**

```bash
# Ensure SISSIz is in PATH
export PATH=$PATH:/path/to/sissiz/bin
which SISSIz
```

**Problem: "R script failed"**

```bash
# Check R packages are installed
R -e 'library(caret); library(randomForest)'

# Verify Rscript location
which Rscript
```

**Problem: Out of memory**

```bash
# Increase Java heap size
java -Xmx16G -jar target/ECSFinder-1.0.0.jar ...
```

**Problem: Very slow processing**

```bash
# Increase CPU cores
-c 16

# Reduce input size or increase thresholds
-mpi 60  # Higher sequence identity threshold
-sszr -4.0  # More stringent structural conservation
```
---

## Repository Structure

```
ECSFinder/
├── src/                          # Java source code
│   └── ca/smithlab/vandalovejoy/ecsfinder/
│       └── ECSFinder.java        # Main class
├── target/
│   └── ECSFinder-1.0.0.jar      # Prebuilt executable JAR
├── RF/                           # Random forest model
│   ├── final_rf_model.rds       # Trained model
│   └── predictions_ECSFinder.R  # Prediction script
├── 3_0_SISSIz/                  # SISSIz resources
├── rfam_families/               # Rfam validation data
├── hybrid_blocks/               # Precomputed alignment blocks
├── example/                     # Example data
│   └── test.maf
├── pom.xml                      # Maven configuration
└── README.md                    # This file
```

---

## Citation

If you use ECSFinder in your research, please cite:
```
Vanda Gaonac’h-Lovejoy, John S Mattick, Martin Sauvageau, Martin A Smith, 
ECSFinder: optimized prediction of evolutionarily conserved RNA secondary structures from genome sequences, 
Nucleic Acids Research, Volume 53, Issue 15, 28 August 2025, gkaf780, https://doi.org/10.1093/nar/gkaf780
```

---

## Acknowledgments

ECSFinder integrates and builds upon several excellent bioinformatics tools and R packages:

#### SISSIz
Gesell T, Washietl S. Dinucleotide controlled null models for comparative RNA gene prediction. *BMC Bioinformatics*. 2008;9:474. doi: 10.1186/1471-2105-9-474

#### ViennaRNA Package (RNALalifold / RNAalifold)
Lorenz R, Bernhart SH, Höner zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL. ViennaRNA Package 2.0. *Algorithms for Molecular Biology*. 2011;6:26. doi: 10.1186/1748-7188-6-26

#### R-scape
Rivas E, Clements J, Eddy SR. A statistical test for conserved RNA structure shows lack of evidence for structure in lncRNAs. *Nucleic Acids Research*. 2017;45(13):8189–8205. doi: 10.1093/nar/gkx386

#### R caret
Kuhn M. Building Predictive Models in R Using the caret Package. *Journal of Statistical Software*. 2008;28(5):1–26. doi: 10.18637/jss.v028.i05

#### R randomForest
Liaw A, Wiener M. Classification and regression by randomForest. *R News*. 2002;2(3):18–22.

---

**Version:** 1.0.0  
**Last Updated:** December 2025
