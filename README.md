


<p align="center" width="100%">
    <img width="25%" src="https://user-images.githubusercontent.com/44384386/195381940-680064be-d53a-45b6-a5e1-a80ff1cb804e.jpg"> 
</p>

# ECSFinder

ECSFinder is a tool designed to scan multiple sequence alignments for **evolutionarily conserved RNA secondary structures**.

Given a set of MAF files, ECSFinder:

1. **Merges and filters alignments** (removing ancestor sequences and duplicates).
2. **Refines local windows** and predicts consensus structures with **RNALalifold**.
3. **Assesses structural conservation** using **SISSIz**, testing whether observed structure is more likely than expected by chance.
4. **Computes structural features** with **RNAalifold** (MFE, pseudo-energy) and sequence conservation metrics (MPI, Shannon entropy, gaps, GC).
5. **Evaluates covariation** using **R-scape**, extracting E-values and base-pair counts for significant helices.
6. **Applies a random forest classifier** (via R) to score each candidate element as likely **true positive (TP)** or **false positive (FP)** based on:
  - E-value
  - Number of significant base pairs
  - Minimal free energy
  - Pseudo-energy
  - Sequence conservation (MPI)
  - Standard deviation of null energies
  - Average MFE from the SISSIz null background
  - Z-score from SISSIz

The result is a robust framework that not only **detects** but also **prioritizes and validates** conserved RNA structures across multi-species alignments. Outputs are suitable for visualization in genome browsers and downstream analysis pipelines.


## Table of Contents

- [Installation](#installation)
  - [SISSIz 3.0](#sissiz-30)
  - [RNALalifold](#rnalalifold)
  - [ECSFinder](#ecsfinder)
  - [R-scape](#r-scape)
  - [R](#r)
- [Usage](#usage)
  - [Running as a class](#running-as-a-class)
  - [Running from JAR](#running-from-jar)
- [Reference species (`-ref`)](#reference-species--ref)
- [Output](#output)
- [Example](#example)


## Installation

### SISSIz 3.0

Install SISSIz according to its README.

Authors:  
Tanja Gesell `<tanja.gesell@univie.ac.at>`  
Stefan Washietl `<wash@tbi.univie.ac.at>`  
Lorenz Perschy `<NA>`

Make sure the `SISSIz` binary is in your `PATH` (ECSFinder calls it via `which SISSIz`).

### RNALalifold

Install the ViennaRNA package from the [ViennaRNA website](https://www.tbi.univie.ac.at/RNA/) and follow the [installation instructions](https://www.tbi.univie.ac.at/RNA/documentation.html#install):

```bash
tar -zxvf ViennaRNA-2.5.1.tar.gz
cd ViennaRNA-2.5.1
./configure
make
sudo make install

Ensure that RNALalifold and RNAalifold are both on your PATH:

### ECSFinder
```
cd ECSFinder/src
javac ECSFinder.java
```
### R-scape

Download the source code [website](http://eddylab.org/R-scape/)
### R

Please use the 4.3 or a more recent version of R. Make sure the caret and randomForest package are installed

## Usage
### ECSFinder
```
java ECSFinder [options] -o output/directory -i input.maf (last parameter must be -i, absolute path required)
 Options:
   -c int number of CPUs for calculations (default 4)
   -g int max gap percentage of sequences for 2D prediction (default 50)
   -sszr double report SISSIzhits below this Z-score (default -3)
   -v verbose (messy but detailed) output
   -mafft realign the maf file before analysis using mafft-ginsi
```

## Output
Three types of results are produced:

* `output.maf` file containing the merged MAF file with the ancestor sequences and duplicate species removed.
* `predicted_ECS.csv` file with predictions made.
* A directory called `aln` containing:
    * A clustal file, e.g., `out_directory/aln/X_9958021_9958096_11_92.2_0.1_0.16_66.8_0.4_304_+.aln`.

### File Name:
***

     1. Name of chromosome (X)
     2. Start loci (9958021)
     3. End loci (9958096)
     4. Number species (11)
     5. Mean pairwise identity (92.2)
     6. Standard deviation (0.1)
     7. Average Shannon entropy (0.16)
     8. GC content (66.8)
     9. Gap content (0.4)
    10. Z-score multiplied by -100 (304)
    11. Direction of the strand (+)
 ***   

*  The genomic coordinates (.bed format) of ECSs are also written to the SDOUT

## Example
 ```
java -jar target/ECSFinder.jar -sszr 0.0 -o TEST -i src/test/resources
chrm	start	end	Num_species	MPI	sd	mean_shannon	gc	gap	zscore	strand	prob
5	RF00017	318	435	5	60.5	0.3	0.73	57.0	1.0	655	+	0.84
5	RF00017	273	426	5	58.1	0.3	0.77	59.2	0.7	1317	+	0.94

...
```