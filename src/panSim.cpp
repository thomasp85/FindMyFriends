#include <Rcpp.h>
#include <iostream>
#include <fstream>

using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix panSim(IntegerVector P, IntegerVector I, CharacterVector names) {
    NumericMatrix res(P.size() - 1, P.size() - 1);
    res.attr("dimnames") = List::create(names, names);
    std::vector<int> setContainer;
    setContainer.reserve(P.size());
    IntegerVector::iterator iterI = I.begin();
    int i;
    int j;
    int k;
    int sim;
    int total;
    int iVal;
    int jVal;
    for(i = 0; i < P.size() - 1; i++) {
        res(i, i) = 1;
        if (i == P.size() - 2) break;
        for(j = i + 1; j < P.size() - 1; j++) {
            setContainer.clear();
            std::set_intersection(iterI + P[i], iterI + P[i + 1] - 1, iterI + P[j], iterI + P[j + 1] - 1, std::back_inserter(setContainer));
            sim = setContainer.size();
            setContainer.clear();
            std::set_union(iterI + P[i], iterI + P[i + 1] - 1, iterI + P[j], iterI + P[j + 1] - 1, std::back_inserter(setContainer));
            total = setContainer.size();
            res(i,j) = double(sim)/total;
            res(j,i) = res(i,j);
        }
    }
    return res;
}
