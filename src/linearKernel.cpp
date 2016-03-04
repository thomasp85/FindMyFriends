#include <Rcpp.h>
#include "fmf-common.h"
using namespace Rcpp;

template<typename T1, typename T2, typename R1, typename R2>
IntegerVector linearKernel(const T1& pX, const T1& jX, const T2& xX, const T1& selX, 
                  double lowerLimit, double upperLimit, R1& P, R1& I, R2& X) {
    int i, j, i1, j1, ind1, ind2, nextMember, nGroups;
    double kv;
    
    int sizeX = selX.size();
    
    // Initialize grouping variables
    std::vector<int> group; // Group representative
    group.reserve(sizeX);
    IntegerVector member(sizeX, -1); // Group membership for each element
    group.push_back(selX[0]); // The first element defines the first group
    I.push_back(0);
    X.push_back(1.0);
    P.push_back(0);
    member[0] = 0;
    std::vector<double> kValues(sizeX, -1); // Temporary storage for kv until it's decided if new group should be created
    
    
    for (j = 1; j < sizeX; j++) {
        R_CheckUserInterrupt();
        
        j1 = selX[j];
        nextMember = -1;
        
        for (i = group.size() - 1; i >= 0; i--) {
            i1 = group[i];
            ind1 = pX[i1];
            ind2 = pX[j1];
            kv = 0;
            
            while (ind1 < pX[i1+1] && ind2 < pX[j1+1]) {
                if (jX[ind1] < jX[ind2]) {
                    ind1++;
                } else if (jX[ind1] > jX[ind2]) {
                    ind2++;
                } else {
                    kv += xX[ind1] * xX[ind2];
                    ind1++;
                    ind2++;
                }
            }
            
            if (kv >= upperLimit) {
                nextMember = i;
                break;
            } else {
                kValues[i] = kv;
            }
        }
        
        if (nextMember != -1) {
            member[j] = nextMember;
        } else {
            P.push_back(I.size());
            nGroups = group.size();
            for(i = 0; i < nGroups; i++) {
                if (kValues[i] > lowerLimit) {
                    I.push_back(i);
                    X.push_back(kValues[i]);
                }
            }
            I.push_back(nGroups);
            X.push_back(1.0);
            group.push_back(j1);
            member[j] = group.size()-1;
        }
    }
    P.push_back(I.size());
    
    return member;
}

//[[Rcpp::export]]
List lkMatrix(IntegerVector pX, IntegerVector jX, NumericVector xX, 
              IntegerVector selX, double lowerLimit, double upperLimit) {
    std::deque<int> P;
    std::deque<int> I;
    std::deque<double> X;
    IntegerVector member = linearKernel(pX, jX, xX, selX, lowerLimit, upperLimit, P, I, X);
    
    return List::create(
        Named("member") = member,
        Named("i") = wrap(I),
        Named("p") = wrap(P),
        Named("x") = wrap(X)
    );
}

//[[Rcpp::export]]
IntegerVector lkMembers(IntegerVector pX, IntegerVector jX, NumericVector xX, 
                        IntegerVector selX, double lowerLimit, 
                        double upperLimit) {
    std::deque<int> P;
    std::deque<int> I;
    std::deque<double> X;
    IntegerVector member = linearKernel(pX, jX, xX, selX, lowerLimit, upperLimit, P, I, X);
    
    IntegerVector grouping = getClusters(I, P, X);
    int i;
    for (i = 0; i < member.size(); ++i) {
        member[i] = grouping[member[i]];
    }
    
    return member;
}
