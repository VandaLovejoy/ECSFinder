


<p align="center" width="100%">
    <img width="25%" src="https://user-images.githubusercontent.com/44384386/195381940-680064be-d53a-45b6-a5e1-a80ff1cb804e.jpg"> 
</p>

ECSFinder is a tool designed to scan multiple alignments for conserved RNA structures. It processes a set of MAF files, calculates key statistics, scans with SISSIz, and outputs BED coordinates of high-confidence predictions. The process begins by refining alignment boundaries using RNALalifold, which identifies locally stable RNA secondary structures. After this refinement, the alignments are analyzed with SISSIz to assess whether a predicted conserved structure is statistically more likely than expected by chance.

Following this, RNAalifold calculates the minimal free energy and pseudo-energy of the predicted structures, providing insights into their stability. R-scape is then used to evaluate the statistical significance of helices within the RNA structures, identifying significant base pairs that further validate the findings.

To enhance prediction accuracy, ECSFinder integrates a random forest model that predicts the likelihood of the identified RNA structures being true positives or false positives. This model considers the following features: E-value, the number of significant base pairs, minimal free energy, pseudo-energy, sequence conservation (MPI), standard deviation, the average MFE from null background from SISSIz and the Z-score.

The result is a robust framework that not only identifies but also validates conserved RNA structures across multiple sequence alignments, providing output that can be visualized and further analyzed using genome browsers and other bioinformatics tools.


## Table of Contents

- [Installation](#installation)
    - [SISSIz](#sissiz)
    - [RNALalifold](#rnalalifold)
    - [ECSFinder](#ecsfinder)
    - [R-scape](#r-scape)
    - [R](#r)
- [Usage](#usage)
- [Output](#output)
- [Example](#example)


## Installation

### SISSIz 3.0

Follow README instructions  
Authors:  
Tanja Gesell <tanja.gesell@univie.ac.at>  
Stefan Washietl <wash@tbi.univie.ac.at>  
Lorenz Perschy <NA>
### RNALalifold
Download the package on the ViennaRNA package [website](https://www.tbi.univie.ac.at/RNA/) and follow the [instructions](https://www.tbi.univie.ac.at/RNA/documentation.html#install)
```
tar -zxvf ViennaRNA-2.5.1.tar.gz
cd ViennaRNA-2.5.1
./configure
make
sudo make install
```

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