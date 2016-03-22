#include <Rcpp.h>
#include "fmf-common.h"

using namespace Rcpp;

//[[Rcpp::export(getClusters)]]
IntegerVector getClustersFromR(IntegerVector I, IntegerVector P, NumericVector X) {
    return getClusters(I, P, X);
}

//[[Rcpp::export]]
List createPanMatrix(IntegerVector org, IntegerVector group) {
    std::deque< std::pair<int, int> > geneRelation;
    std::deque<int> I, P, X;
    int i;
    int currentOrg = 0;
    int currentGroup = 0;
    
    for (i = 0; i < org.size(); ++i) {
        geneRelation.push_back(std::pair<int, int>(org[i], group[i]));
    }
    std::sort(geneRelation.begin(), geneRelation.end());
    
    for (i = 0; i < geneRelation.size(); ++i) {
        while (currentOrg != geneRelation[i].first) {
            P.push_back(I.size());
            currentGroup = 0;
            ++currentOrg;
        }
        if (currentGroup == geneRelation[i].second) {
            ++X.back();
        } else {
            currentGroup = geneRelation[i].second;
            I.push_back(currentGroup);
            X.push_back(1);
        }
    }
    P.push_back(I.size());
    
    return List::create(
        Named("i") = wrap(I),
        Named("p") = wrap(P),
        Named("x") = wrap(X)
    );
}

//[[Rcpp::export]]
DataFrame calcGroupInfo(List groupOrgs, int nOrgs, double threshold) {
    IntegerVector geneToOrg, uGeneToOrg;
    CharacterVector group(groupOrgs.size());
    IntegerVector nOrg(groupOrgs.size());
    IntegerVector nGenes(groupOrgs.size());
    
    for (int i = 0; i < groupOrgs.size(); ++i) {
        geneToOrg = groupOrgs[i];
        uGeneToOrg = unique(geneToOrg);
        if (uGeneToOrg.size() / double(nOrgs) >= threshold) {
            group[i] = "Core";
        } else if (uGeneToOrg.size() == 1) {
            group[i] = "Singleton";
        } else {
            group[i] = "Accessory";
        }
        nOrg[i] = uGeneToOrg.size();
        nGenes[i] = geneToOrg.size();
    }
    
    return DataFrame::create(
        Named("group") = group,
        Named("nOrg") = nOrg,
        Named("nGenes") = nGenes,
        Named("stringsAsFactors") = false
    );
}

//[[Rcpp::export]]
IntegerVector findIn(IntegerVector keys, IntegerVector lookup) {
    std::deque<int> res;
    int i, j;
    int nl = lookup.size();
    int nk = keys.size();
    for (i = 0; i < nl; ++i) {
        for (j = 0; j < nk; ++j) {
            if (lookup[i] == keys[j]) {
                res.push_back(i+1);
                break;
            }
        }
    }
    return wrap(res);
}
