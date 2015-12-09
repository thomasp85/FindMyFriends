#include <Rcpp.h>
using namespace Rcpp;

// static int *ptrP;
// static int *ptrI;
// static double *ptrX;

//[[Rcpp::export]]
List linearKernel(IntegerVector pX, IntegerVector jX, 
                  NumericVector xX, IntegerVector selX,
                  double lowerLimit, double upperLimit) {
    int i, j, i1, j1, ind1, ind2, nextMember, nGroups;
    double kv;
    
    int sizeX = selX.size();
    
    // set initial size and growth factor for i and x arrays
    uint64_t vectorSize = sizeX * 2;
//     double growBy = 1.6;
//     
//     // int *pptr = (int *) Calloc(sizeX + 1, int);
//     int *iptr = (int *) Calloc(vectorSize, int);
//     double *xptr = (double *) Calloc(vectorSize, double);
//     // ptrP = pptr;
//     ptrI = iptr;
//     ptrX = xptr;
//     
//     int nextFree = 0;
    
    std::deque<int> P;
    std::deque<int> I;
    std::deque<double> X;
//     P.reserve(sizeX + 1);
//     I.reserve(vectorSize);
//     X.reserve(vectorSize);
    
    // Initialize grouping variables
    std::vector<int> group; // Group representative
    group.reserve(sizeX);
    IntegerVector member(sizeX, -1); // Group membership for each element
    group.push_back(selX[0]); // The first element defines the first group
    I.push_back(0);
    X.push_back(1.0);
    P.push_back(0);
    // pptr[nextFree] = 0;
//     iptr[nextFree] = 0;
//     xptr[nextFree++] = 0;
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
            // P.push_back(nextFree);
            nGroups = group.size();
//             if (((uint64_t) nextFree + nGroups) > vectorSize) {
//                 // realloc arrays i and x
//                 vectorSize += nGroups;
//                 vectorSize *= growBy;
//                 iptr = (int *) Realloc(iptr, vectorSize, int);
//                 xptr = (double *) Realloc(xptr, vectorSize, double);
//                 ptrI = iptr;
//                 ptrX = xptr;
//             }
            for(i = 0; i < nGroups; i++) {
                if (kValues[i] > lowerLimit) {
                    I.push_back(i);
                    X.push_back(kValues[i]);
//                     iptr[nextFree] = i;
//                     xptr[nextFree++] = kValues[i];
                }
            }
//             iptr[nextFree] = nGroups;
//             xptr[nextFree++] = 1.0;
            I.push_back(nGroups);
            X.push_back(1.0);
            group.push_back(j1);
            member[j] = group.size()-1;
        }
    }
    P.push_back(I.size());
//     P.push_back(nextFree);
//     
//     
//     IntegerVector I;
//     NumericVector X;
//     if (nextFree > 0) {
//         I = IntegerVector(iptr, iptr + nextFree);
//         X = NumericVector(xptr, xptr + nextFree);
//     }
//     
//     // Release memory
//     Free(iptr);
//     ptrI = NULL;
//     Free(xptr);
//     ptrX = NULL;
    
    return List::create(
        Named("member") = member,
        Named("i") = wrap(I),
        Named("p") = wrap(P),
        Named("x") = wrap(X)
    );
}

// //[[Rcpp::export]]
// void freeHeapLinearKernelC() {
//     if (ptrP != NULL)
//     {
//         Free(ptrP);
//         ptrP = NULL;
//     }
//     
//     if (ptrI != NULL)
//     {
//         Free(ptrI);
//         ptrI = NULL;
//     }
//     
//     if (ptrX != NULL)
//     {
//         Free(ptrX);
//         ptrX = NULL;
//     }
// }

//[[Rcpp::export]]
    
    }
    
}
