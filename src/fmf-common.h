#ifndef FMF_COMMON
#define FMF_COMMON

#include <Rcpp.h>

using namespace Rcpp;

template<typename T1, typename T2>
IntegerVector getClusters(const T1& I, const T1& P, const T2& X) {
    int currentCluster = 1;
    int node;
    int i;
    int nextNode;
    int nNodes = P.size() - 1;
    IntegerVector members(nNodes, 0);
    std::vector<bool> visited(nNodes, false);
    std::vector<int> nodeStack;
    nodeStack.reserve(nNodes);
    int seed = 0; // Set first seed
    nodeStack.push_back(seed);
    members[seed] = currentCluster;
    
    while(true) {
        do {
            node = nodeStack.back();
            nodeStack.pop_back();
            if (visited[node]) {
                continue;
            } else {
                visited[node] = true;
            }
            i = P[node];
            while (i < P[node+1]) {
                nextNode = I[i++];
                if (members[nextNode] == 0) {
                    members[nextNode] = currentCluster;
                    nodeStack.push_back(nextNode);
                }
            }
        } while(nodeStack.size());
        
        while (members[seed]) {
            if (++seed >= nNodes) {
                break;
            }
        }
        if (seed >= nNodes) break;
        
        nodeStack.push_back(seed);
        members[seed] = ++currentCluster;
    }
    
    return members;
}

#endif
