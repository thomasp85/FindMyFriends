#include <Rcpp.h>
#include "fmf-common.h"

using namespace Rcpp;

//[[Rcpp::export(getClusters)]]
IntegerVector getClustersFromR(IntegerVector I, IntegerVector P, NumericVector X) {
    return getClusters(I, P, X);
}
