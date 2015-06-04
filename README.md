# FindMyFriends
*Fast alignment-free pangenome creation and exploration*

[![Build Status](https://travis-ci.org/thomasp85/FindMyFriends.svg?branch=master)](https://travis-ci.org/thomasp85/FindMyFriends)

FindMyFriends is an R package for doing pangenomic analyses on microbial genomes. It is slated for inclusion in Bioconductor but until then it can be installed using devtools:

```R
if(!require(devtools)) {
  install.packages('devtools')
  library(devtools)
}
install_github('thomasp85/FindMyFriends')
```

## Methodology
A lot of different tools exists for calculating pangenomes (PanOCT, OrthoMCL, PanFunPro, GET_HOMOLOGUES etc.). Most of these (PanFunPro excepted) utilize BLAST as a core step for generating the data needed for grouping genes. BLAST is great and has had a profound effect on the world of biology and bioinformatics, but it is slow. Consequently generation of pangenomes using BLAST based tools quickly becomes unfeasible on standard hardware as the number of genes in the pangenome increases. FindMyFriends circumvent the time consuming BLAST analysis be relying on kmer-derived cosine similarities. Furthermore it utilizes a new heuristic approach to avoid comparing all genes with each other, thus avoiding the quadratic scaling problem that has riddled pangenome analysis for years.
