#include <Rcpp.h>
#include <iostream>
#include <fstream>

using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix panSim(NumericMatrix pg) {
    NumericMatrix res(pg.ncol(), pg.ncol());
    List dimnames = pg.attr("dimnames");
    res.attr("dimnames") = List::create(dimnames[1], dimnames[1]);
    int i;
    int j;
    int k;
    int sim;
    int total;
    int iVal;
    int jVal;
    for(i=0; i < pg.ncol(); i++) {
        res(i,i) = 1;
        for(j=i+1; j < pg.ncol(); j++) {
            sim = 0;
            total = 0;
            for(k=0; k < pg.nrow(); k++) {
                iVal = pg(k,i);
                jVal = pg(k,j);
                if(iVal == 0 && jVal == 0) continue;
                if(iVal != 0 && jVal != 0) sim++;
                total++;
            }
            res(i,j) = double(sim)/total;
            res(j,i) = res(i,j);
        }
    }
    return res;
}