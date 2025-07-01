# Rfam Subset Generator & Hybrid-Block Alignments  

This repository contains the sequence subsets and hybrid MAF blocks that were
used to train and test our RNA-structure classifiers.  
All datasets were produced *de-novo* with the **GenerateRFAMsubsets** tool

Check out the original project:
> [noncodo/BigRedButton](https://github.com/noncodo/BigRedButton).

---

## 1 . Overview

``GenerateRFAMsubsets`` walks through an input folder full of **_aligned_**
FASTA files from [Rfam](https://rfam.org) and picks a *balanced* sample of
sequences per family.

---

## 2 . Building the tool

```bash
cd path/to/GenerateRFAMsubsets
javac com/example/oldrfam/GenerateRFAMsubsets.java