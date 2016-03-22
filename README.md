# FindMyFriends
*Fast alignment-free pangenome creation and exploration*

[![How long since the package was first in a released Bioconductor version (or is it in devel only).](http://bioconductor.org/shields/years-in-bioc/FindMyFriends.svg)](http://bioconductor.org/packages/FindMyFriends)
[![Percentile (top 5/20/50% or 'available') of downloads over last 6 full months. Comparison is done across all package categories (software, annotation, experiment).](http://bioconductor.org/shields/downloads/FindMyFriends.svg)](http://bioconductor.org/packages/stats/bioc/FindMyFriends.html) Release: 
[![Build status in Release branch.](http://bioconductor.org/shields/build/release/bioc/FindMyFriends.svg)](http://bioconductor.org/packages/release/bioc/html/FindMyFriends.html) Devel: 
[![Build status in Devel branch.](http://bioconductor.org/shields/build/devel/bioc/FindMyFriends.svg)](http://bioconductor.org/packages/devel/bioc/html/FindMyFriends.html) @master: 
[![Build Status](https://travis-ci.org/thomasp85/FindMyFriends.svg?branch=master)](https://travis-ci.org/thomasp85/FindMyFriends) 
[![codecov.io](http://codecov.io/github/thomasp85/FindMyFriends/coverage.svg?branch=master)](http://codecov.io/github/thomasp85/FindMyFriends?branch=master)

FindMyFriends is an R package for doing pangenomic analyses on microbial 
genomes. It is released as part of the [Bioconductor](http://bioconductor.org/) 
project and can be installed with the `biocLite()` function:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("FindMyFriends")
```

For the absolute latest version, install directly from GitHub:

```R
if(!require(devtools)) {
  install.packages('devtools')
  library(devtools)
}
install_github('thomasp85/FindMyFriends')
```

## But what is this really about?
In comparative microbial genomics a pangenome is defined as a grouping of genes
across genomes based on some sort of similarity. This similarity meassure is not
set in stone, but often it is derived from BLASTing each pair of genes against 
each other. This is a bad idea for several reasons: comparing all against all
leads to a horrible scaling of computational time as the number of genes 
increase, BLAST is in general really slow, and sequence similarity alone cannot
distinguish orthologue genes from paralogues. The last point has been adressed
by recent tools such as PanOCT and Roary, but the first two still stands (though
Roary do something clever to make it less of an issue).

Enter FindMyFriends...

## So this is just another algorithm?
It is also another algorithm. But more importantly it is a framework for 
conducting pangenome analysis that is completely agnostic to how you've derived
your pangenome in the first place. FindMyFriends defines an extendible list of
classes for handling pangenome data in a transparent way, and plugs directly 
into the vast array of genomic tools offered by Bioconductor.

Okay, back to the algorithms. FindMyFriends works by using CD-Hit to create a 
very coarse grouping of the genes in your dataset, and then refine this grouping
in a second pass using additional similarity meassures. This is in contrast to
Roary that uses CD-Hit, but only to group the most similar genes together prior
to running BLAST. The second pass in FindMyFriends is where all the magic 
happens. The genes in each large group is compared by sequence similarity (using
kmer cosine similarity), sequence length, genome membership and neighborhood 
similarity. Based on these comparisons a graph is created for each group, with 
edges defining similarity above a certain threshold between genes. From this 
graph cliques are gradually extracted in a way that ensures the highest quality
cliques are extracted first. These cliques defines the final grouping of genes.
Because they are cliques the user can be sure that all members of the resulting
gene groups share a defined similarity with each other and that no gene can be
grouped with others solemnly based on a high similarity to one member.

## That sound like a lot of work
Well, high quality results are more important than speed! But this is one of the
rare cases where you can have your cake and eat it too. FindMyFriends is, by a 
large margin, the fastest algorithm out there:

![Comparison of computational time for different pangenome algorithms.](https://dl.dropboxusercontent.com/u/2323585/FindMyFriends/timing.png)

FindMyFriends scales to thousands of genomes, and can handle large diversity 
(i.e. not restricted to species level). As an example a pangenome based on ~1200 
strains from the order Lactobacillales (Lactic Acid Bacteria) was created in 
around 8 hours on a c3x8.large AWS instance using a single core.

## How do I use it then?
Being a framework there is a lot of things you can do and many different ways
to do it. Following is the recommended approach to calculating a pangenome:

```r
library(FindMyFriends)

# We expect here that your genomes are stored in amino acid fasta files in the
# working directory.
genomes <- list.files(pattern = '.fasta')

# First we create our pangenome object
pg <- pangenome(genomes, translated = TRUE, geneLocation = 'prodigal')

# Then we make the initial grouping
pg <- cdhitGrouping(pg)

# And lastly we refine the groups
pg <- neighborhoodSplit(pg)
```

please see the vignette for more information on the different steps as well as
examples on what you can do with your data once you're done grouping your genes.

## What will happen with this in the future?
Following are some of the features that are being worked on/considered:

- GFF3 and GBK file support
- Improved plotStat and plotEvolution
- Even more panchromosomal analysis tools
    1. plotPC
    2. Automatic frameshift detection
    3. Support for storing pc-derived grouping in object
- Better parallelization
- Switch to using data.table internally for better performance
- Exporting functions
- sqlite based class for very low memory interface
