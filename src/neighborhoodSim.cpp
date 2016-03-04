

#include <Rcpp.h>
#include <algorithm>
#include "fmf-common.h"
#include "progress.h"

using namespace Rcpp;

void getNeighbors(IntegerVector up, IntegerVector down, IntegerVector groups, int start, int length, std::vector<int> & result) {
    std::vector<int> neighborsUp;
    std::vector<int> neighborsDown;
    int next, i;
    
    next = up[start];
    for (i = 0; i < length; i++) {
        if (next == -1) break;
        
        neighborsUp.push_back(groups[next]);
        next = up[next];
    }
    next = down[start];
    for (i = 0; i < length; i++) {
        if (next == 0) break;
        
        neighborsDown.push_back(groups[next]);
        next = down[next];
    }
    result.clear();
    result.insert(result.end(), neighborsDown.rbegin(), neighborsDown.rend());
    result.insert(result.end(), neighborsUp.begin(), neighborsUp.end());
}

void filterNeighborhood(std::vector<int> neighborhood, std::vector<int> lookup, std::vector<int> & result) {
    for (int i = 0; i < neighborhood.size(); i++) {
        if (std::find(lookup.begin(), lookup.end(), neighborhood[i]) != lookup.end())
            result.push_back(neighborhood[i]);
    }
}

//[[Rcpp::export]]
List neighborhoodSim(IntegerVector members, IntegerVector groups, 
                     IntegerVector organism, int size, IntegerVector down, 
                     IntegerVector up, LogicalVector reverse, 
                     IntegerVector width, double threshold, 
                     bool forceParalogues) {
    int i, j, k, i1, j1, widthi, widthj, sim;
    double diff;
    
    int nMembers = members.size();
    
    // Neighborhood storage
    std::vector< std::vector<int> > neighborhood(nMembers, std::vector<int>(size*2));
    std::vector<bool> hasNeighbors(nMembers, false);
    std::vector<int> neighborhoodI, neighborhoodJ;
    
    std::vector<int> P;
    std::deque<int> I;
    std::deque<int> X;
    P.reserve(nMembers + 1);
    
    for (j = 0; j < nMembers - 1; j++) {
        // R_CheckUserInterrupt();
        
        j1 = members[j];
        
        P.push_back(I.size());
        
        for (i = j + 1; i < nMembers; i++) {
            i1 = members[i];
            
            // Check paralogue
            if (forceParalogues && organism[i1] == organism[j1])
                continue;
            
            // Check sequence lengths
            widthi = width[i1];
            widthj = width[j1];
            diff = abs(widthi - widthj);
            if (threshold < 1) {
                if (diff/std::min(widthi, widthj) > threshold)
                    continue;
            } else {
                if (diff > threshold)
                    continue;
            }
            
            // Calculate similarity
            neighborhoodI.clear();
            neighborhoodJ.clear();
            sim = 0;
            
            if (!hasNeighbors[i]) {
                if (reverse[i1]) {
                    getNeighbors(down, up, groups, i1, size, neighborhood[i]);
                } else {
                    getNeighbors(up, down, groups, i1, size, neighborhood[i]);
                }
                hasNeighbors[i] = true;
            }
            if (!hasNeighbors[j]) {
                if (reverse[j1]) {
                    getNeighbors(down, up, groups, j1, size, neighborhood[j]);
                } else {
                    getNeighbors(up, down, groups, j1, size, neighborhood[j]);
                }
                hasNeighbors[j] = true;
            }
            filterNeighborhood(neighborhood[i], neighborhood[j], neighborhoodI);
            filterNeighborhood(neighborhood[j], neighborhoodI, neighborhoodJ);
            
            for (k = 0; k < neighborhoodI.size(); k++) {
                if (neighborhoodI[k] == neighborhoodJ[k])
                    sim++;
            }
            if (sim != 0) {
                I.push_back(i);
                X.push_back(sim);
            }
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
DataFrame mergeSims(IntegerVector nI, IntegerVector nP, IntegerVector nX, 
                    IntegerVector sI, IntegerVector sP, NumericVector sX, 
                    IntegerVector guideGroup) {
    int maxSize = nX.size() > sX.size() ? nX.size() : sX.size();
    std::vector<int> from, to, nSim, gSim;
    std::vector<double> sSim;
    from.reserve(maxSize);
    to.reserve(maxSize);
    nSim.reserve(maxSize);
    gSim.reserve(maxSize);
    sSim.reserve(maxSize);
    int i, i1, i2, j1, j2;
    
    for (i = 0; i < nP.size() - 1; i++) {
        R_CheckUserInterrupt();
        
        i1 = nP[i];
        i2 = sP[i];
        while (i1 < nP[i+1] && i2 < sP[i+1]) {
            j1 = nI[i1];
            j2 = sI[i2];
            if (j1 > j2) {
                i2++;
            } else if (j1 < j2) {
                i1++;
            } else {
                from.push_back(i+1);
                to.push_back(j1+1);
                nSim.push_back(nX[i1]);
                sSim.push_back(sX[i2]);
                gSim.push_back(guideGroup[i] == guideGroup[j1] ? 1 : 0);
                i1++;
                j1++;
            }
        }
    }
    
    return DataFrame::create(
        Named("from") = wrap(from),
        Named("to") = wrap(to),
        Named("nSim") = wrap(nSim),
        Named("sSim") = wrap(sSim),
        Named("gSim") = wrap(gSim)
    );
}

// //[[Rcpp::export]]
// IntegerVector widthSim(IntegerVector groups, IntegerVector width, double threshold) {
//     int i, j, minLength, widthi, widthj;
//     double diff;
//     
//     int size = groups.size();
//     
//     std::deque<int> P;
//     std::deque<int> I;
//     std::deque<int> X;
//     
//     for (j = 0; j < size - 1; j++) {
//         R_CheckUserInterrupt();
//         
//         P.push_back(I.size());
//         for (i = j+1; i < size; i++) {
//             if (groups[j] != groups[i])
//                 continue;
//             
//             // Check sequence lengths
//             widthi = width[i];
//             widthj = width[j];
//             minLength = widthi < widthj ? widthi : widthj;
//             diff = abs(widthi - widthj);
//             if (threshold < 1) {
//                 if (diff/minLength > threshold)
//                     continue;
//             } else {
//                 if (diff > threshold)
//                     continue;
//             }
//             I.push_back(i);
//             X.push_back(1);
//         }
//     }
//     P.push_back(I.size());
//     P.push_back(I.size());
//     
//     IntegerVector newgroups = getClusters(I, P, X);
//     return newgroups;
// }

//[[Rcpp::export]]
IntegerVector widthSim(List groups, IntegerVector width, double threshold, CharacterVector progName, bool showProgress) {
    IntegerVector res(width.size());
    int i, j, k, nMembers, id1, id2, widthi, widthj;
    int maxgroup = 0;
    double diff;
    
    int size = groups.size();
    
    std::deque<int> P;
    std::deque<int> I;
    std::deque<int> X;
    
    Progress prog(size + 1, as<std::string>(progName), 100, showProgress);
    prog.start();
    
    for (k = 0; k < size; ++k) {
        IntegerVector group = groups[k];
        P.clear();
        I.clear();
        X.clear();
        nMembers = group.size();
        if (nMembers == 1) {
            res[group[0] - 1] = maxgroup + 1;
            ++maxgroup;
        } else {
            for (j = 0; j < nMembers - 1; ++j) {
                id1 = group[j] - 1;

                P.push_back(I.size());
                for (i = j+1; i < nMembers; ++i) {
                    id2 = group[i] - 1;

                    // Check sequence lengths
                    widthi = width[id2];
                    widthj = width[id1];
                    diff = abs(widthi - widthj);
                    if (threshold < 1) {
                        if (diff/std::min(widthi, widthj) > threshold)
                            continue;
                    } else {
                        if (diff > threshold)
                            continue;
                    }
                    I.push_back(i);
                    X.push_back(1);
                }
            }
            P.push_back(I.size());
            P.push_back(I.size());

            IntegerVector newgroups = getClusters(I, P, X);
            newgroups = newgroups + maxgroup;
            maxgroup = max(newgroups);
            for (i = 0; i < nMembers; ++i) {
                res[group[i] - 1] = newgroups[i];
            }
        }
        prog.increment();
    }
    prog.finish();
    
    return res;
}

//[[Rcpp::export]]
IntegerVector getCliques(RObject graph) {
    // R functionality
    Environment FMF("package:FindMyFriends");
    List igraph = FMF[".igraphFunctions"];
    Function neighbors = igraph["neighbors"];
    Function gorder = igraph["gorder"];
    Function gsize = igraph["gsize"];
    Function ends = igraph["ends"];
    Function vertexAttr = igraph["vertex_attr"];
    Function deleteVertices = igraph["delete_vertices"];
    
    // Containers for graph information etc.
    int nVertices = as<int>(gorder(graph));
    int nEdges = as<int>(gsize(graph));
    IntegerMatrix edges;
    int i;
    int cliqueID = 1;
    
    // Result variable
    IntegerVector cliques(nVertices);
    
    // Storage of possible vertices to add
    int * possibles = new int[nVertices];
    int * possiblesNext = new int[nVertices];
    int * possiblesTemp;
    int * possiblesEnd;
    
    // Current clique storage and touched edge info
    std::set<int> cliqueMembers;
    std::vector<bool> disregard;
    disregard.reserve(nEdges);
    
    // Initialize possibles
    IntegerVector neighbor1, neighbor2;
    
    while (nEdges) {
        R_CheckUserInterrupt();
        
        edges = ends(graph, seq_len(nEdges), wrap(false));
        IntegerMatrix::Column from = edges(_, 0);
        IntegerMatrix::Column to = edges(_, 1);
        
        cliqueMembers.insert(from[0]);
        cliqueMembers.insert(to[0]);
        
        neighbor1 = neighbors(graph, from[0]);
        neighbor2 = neighbors(graph, to[0]);
        possiblesEnd = std::set_intersection(neighbor1.begin(), neighbor1.end(), neighbor2.begin(), neighbor2.end(), possibles);
        
        disregard.resize(nEdges);
        std::fill(disregard.begin(), disregard.end(), false);
        
        while (possibles != possiblesEnd) {
            for (i = 1; i <= nEdges; i++) {
                if (i == nEdges) // Couldn't find a qualifying edge while possibles remain
                    stop("Missing edge to possible member");
                if (disregard[i])
                    continue;
                if (std::binary_search(possibles, possiblesEnd, from[i])) {
                    if (cliqueMembers.find(to[i]) == cliqueMembers.end())
                        continue;
                    neighbor1 = neighbors(graph, from[i]);
                    cliqueMembers.insert(from[i]);
                } else if (std::binary_search(possibles, possiblesEnd, to[i])) {
                    if (cliqueMembers.find(from[i]) == cliqueMembers.end())
                        continue;
                    neighbor1 = neighbors(graph, to[i]);
                    cliqueMembers.insert(to[i]);
                } else {
                    disregard[i] = true;
                    continue;
                }
                disregard[i] = true;
                possiblesEnd = std::set_intersection(possibles, possiblesEnd, neighbor1.begin(), neighbor1.end(), possiblesNext);
                
                // Swap possibles storage around
                possiblesTemp = possibles;
                possibles = possiblesNext;
                possiblesNext = possiblesTemp;
                
                break;
            }
        }
        IntegerVector members(cliqueMembers.begin(), cliqueMembers.end());
        IntegerVector memberIDs = vertexAttr(graph, "ID", members);
        for (IntegerVector::iterator memberIt = memberIDs.begin(); memberIt != memberIDs.end(); ++memberIt) {
            cliques[*memberIt] = cliqueID;
        }
        cliqueID++;
        cliqueMembers.clear();
        graph = deleteVertices(graph, members);
        nEdges = as<int>(gsize(graph));
    }
    
    if (gorder(graph) != 0) {
        for (i = 0; i < cliques.size(); ++i) {
            if (cliques[i] == 0) {
                cliques[i] = cliqueID;
                cliqueID++;
            }
        }
    }
    delete[] possibles;
    delete[] possiblesNext;
    return cliques;
}

//[[Rcpp::export]]
IntegerVector getPotentials(IntegerVector down, IntegerVector up, 
                            LogicalVector pending, LogicalVector reverse,
                            List groupSplit, IntegerVector groups) {
    std::vector<int> potentials;
    potentials.reserve(pending.size());
    IntegerVector group;
    int i, j, j1, previous, next;
    bool before, after;
    
    for (i = 0; i < pending.size(); ++i) {
        if (pending[i]) {
            group = groupSplit[i];
            before = after = true;
            for (j = 0; j < group.size(); ++j) {
                if (!after && !before)
                    break;
                j1 = group[j]-1;
                if (before) {
                    previous = reverse[j1] ? up[j1] : down[j1];
                    if (previous != -1 && pending[groups[previous]]) {
                        before = false;
                    }
                }
                if (after) {
                    next = reverse[j1] ? down[j1] : up[group[j]-1];
                    if (next != -1 && pending[groups[next]]) {
                        after = false;
                    }
                }
            }
            if (before || after) {
                potentials.push_back(i + 1);
            }
        }
    }
    
    return wrap(potentials);
}

//[[Rcpp::export]]
LogicalVector groupHasParalogues(List groupMembers, IntegerVector org) {
    int nGroups = groupMembers.size();
    int i, j, currentOrg;
    int maxO = 0;
    IntegerVector group;
    LogicalVector hasParalogues(nGroups, false);
    std::vector<bool> visited(max(org), false);
    
    for (i = 0; i < nGroups; ++i) {
        R_CheckUserInterrupt();
        
        group = groupMembers[i];
        for (j = 0; j < group.size(); ++j) {
            currentOrg = org[group[j]-1]-1;
            if (currentOrg < visited.size() && visited[currentOrg]) {
                hasParalogues[i] = true;
                break;
            } else {
                visited[currentOrg] = true;
            }
        }
        while (j != 0) {
            visited[org[group[--j]-1]-1] = false;
        }
    }
    
    return hasParalogues;
}
