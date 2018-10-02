# FastaUtils: Utilities for DNA/RNA sequence processing


The package **FastaUtils** provides tools for the manipulation of DNA/RNA sequence data with R. This package is still *work in progress*.

**FastaUtils** depends on **Biostrings** package which can be installed from **Bioconductor** (see below).

## Content
The current version contains the following functions:

+ **subset.fasta()** 
Select a subset of sequences from a fasta file
+  **DNA2RNA()** 
Convert a DNA fasta to an RNA fasta
+ **RNA2DNA()**
Convert a RNA fasta to an DNA fasta
+ **fasta.sample()**
Select a random sample of sequences
+ **fasta.cutter()**
Cuts peaces of a fixed length and a random position for all sequences or for a random sample of sequences

## Authors

Guillem Salazar

## Citation

To see the preferable citation of the package, type:

```r
citation("FastalUtils")
```
## Software Versions

*R*: version 3.2.1  
*RStudio*: Version 0.99.467

## Contact Info

Guillem Salazar (guillems at ethz.ch)

### Installation of the package

* Install the latest version of **Biostring** from Bioconductor:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
```

* Install FastaUtils' current development version from Github:

```r
devtools::install_github("GuillemSalazar/FastaUtils")
```

* After attaching the package, you are ready to start:

```r
library(FastaUtils)
```


