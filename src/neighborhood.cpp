//#include <algorithm>
//#include <set>
//
//using namespace Rcpp;
//using namespace std;
//
////[[Rcpp::export]]
//NumericMatrix neighborhoodSim(List down, List up) {
//    if(length(down) != length(up))
//        stop("Length of both input list must be equal");
//        
//    NumericMatrix dist(length(down), length(down));
//    for(int i=0; i<length(down); i++) {
//        for(int j=i+1; j<length(down); j++) {
//            set<string> downSet = myIntersect(down[i], down[j]);
//            if(downSet.empty()) continue;
//            
//            set<string> upSet = myIntersect(up[i], up[j]);
//            if(upSet.empty()) continue;
//        }
//    }
//    return(dist);
//}
//
//set<string> myIntersect(CharacterVector a, CharacterVector b) {
//    set<string> setA(a.begin(), a.end());
//    set<string> setB(b.begin(), b.end());
//    set<string> result;
//    set_intersection(setA.begin(), setA.end(), setB.begin(), setB.end(), inserter(result,result.begin()));
//    return(result);
//}