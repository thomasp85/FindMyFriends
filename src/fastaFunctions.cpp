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
    }
    return counter;
}
struct namedSequences {
    std::vector<std::string> names;
    std::vector<std::string> sequences;
};
namedSequences readSequences(std::string filename, std::vector<int> index) {
    const int BUFFERSIZE = 5000;
    const int NSEQS_BUFFER = 10000;
    std::vector<std::string> sequences;
    std::vector<std::string> names;
    bool all = false;
    if(index[0] != -1) {
        sequences.reserve(index.size());
        names.reserve(index.size());
    } else {
        sequences.reserve(NSEQS_BUFFER);
        names.reserve(NSEQS_BUFFER);
        all = true;
    }
    int counter = 0;
    bool reading = false;
    std::string currentSequence = "";
    currentSequence.reserve(BUFFERSIZE);
    std::vector<int>::iterator locator;
    std::string line;
    std::ifstream in(filename.c_str());
    if(in.is_open()) {
        while(getline(in,line)) {
            if (!line.empty() && line[line.size() - 1] == '\r') {
                line.erase(line.size() - 1);
            }
            if(line[0] == '>') {
                counter++;
                if(reading) {
                    sequences.push_back(currentSequence);
                    currentSequence = "";
                    reading = false;
                    if(!all) {
                        if(index.size() == 0) break;
                    }
                }
                locator = std::find(index.begin(), index.end(), counter);
                if(all || locator != index.end()) {
                    names.push_back(line.substr(1));
                    reading = true;
                    if(!all) {
                        index.erase(locator);
                    }
                    continue;
                }
            }
            if(reading) {
                currentSequence += line;
            }
        }
        if(reading) {
            sequences.push_back(currentSequence);
        }
    }
    namedSequences res = {
        names,
        sequences
    };
    return res;
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

//[[Rcpp::export]]
List readFasta(CharacterVector files, List indexes) {
    int nFiles = files.size();
    List res(nFiles);
    namedSequences tempRes;
    
    for(int i=0; i<nFiles; i++) {
        tempRes = readSequences(as<std::string>(files[i]), as< std::vector<int> >(indexes[i]));
        CharacterVector tempSeq = wrap(tempRes.sequences);
        CharacterVector tempNames = wrap(tempRes.names);
        tempSeq.names() = tempNames;
        res[i] = tempSeq;
    }
    
    return res;
}

//[[Rcpp::export]]
List test() {
    List res(2);
    return res;
}
