#include <Rcpp.h>
#include <iostream>
#include <fstream>

using namespace Rcpp;

// Count number of '>' as a proxy for the number of sequences in a fasta file
int countSeq(std::string filename) {
    int counter = 0;
    std::string line;
    std::ifstream in(filename.c_str());
    if(in.is_open()) {
        while(getline(in,line)) {
            if(line[0] == '>') counter++;
        }
        in.close();
    }
    return counter;
}

//[[Rcpp::export]]
NumericVector nSeqs(CharacterVector files) {
    int nFiles = files.size();
    NumericVector res(nFiles);
    
    for(int i=0; i<nFiles; i++) {
        res[i] = countSeq(as<std::string>(files[i]));
    }
    
    return res;
}

