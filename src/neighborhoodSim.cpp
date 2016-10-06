

#include <Rcpp.h>
#include <algorithm>
#include "fmf-common.h"
#include "smallgraph.h"
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
IntegerVector getCliques(DataFrame edges, int nNodes) {
    Graph gr(nNodes, edges);
    
    if (gr.isComplete()) {
        return IntegerVector(nNodes, 1);
    }
    
    // Containers for graph information etc.
    int nEdges = edges.nrows();
    int i;
    int cliqueID = 1;
    
    // Result variable
    IntegerVector cliques(nNodes);
    
    // Storage of possible vertices to add
    std::vector<int> possibles, possiblesNext;
    possibles.reserve(nNodes);
    possiblesNext.reserve(nNodes);
    
    // Current clique storage and touched edge info
    std::vector< std::pair<int, int> > edgelist;
    std::vector<int> neighbor1, neighbor2, members;
    std::vector<int>::iterator itvec;
    std::set<int> cliqueMembers;
    std::vector<bool> disregard;
    disregard.reserve(nEdges);
    
    while (nEdges) {
        R_CheckUserInterrupt();
        
        if (gr.isComplete()) {
            members = gr.nodeIds();
        } else {
            edgelist = gr.completeTriangle(gr.firstEdge());
            
            cliqueMembers.insert(edgelist[0].first);
            cliqueMembers.insert(edgelist[0].second);
            
            if (edgelist.size() > 1) {
                neighbor1 = gr.neighbors(edgelist[0].first);
                neighbor2 = gr.neighbors(edgelist[0].second);
                possibles.clear();
                std::set_intersection(neighbor1.begin(), neighbor1.end(), neighbor2.begin(), neighbor2.end(), std::back_inserter(possibles));
                
                disregard.resize(nEdges);
                std::fill(disregard.begin(), disregard.end(), false);
                
                while (possibles.size() != 0) {
                    for (i = 1; i <= nEdges; i++) {
                        if (i == nEdges) // Couldn't find a qualifying edge while possibles remain
                            stop("Missing edge to possible member");
                        if (disregard[i])
                            continue;
                        if (std::binary_search(possibles.begin(), possibles.end(), edgelist[i].first)) {
                            if (cliqueMembers.find(edgelist[i].second) == cliqueMembers.end())
                                continue;
                            neighbor1 = gr.neighbors(edgelist[i].first);
                            cliqueMembers.insert(edgelist[i].first);
                        } else if (std::binary_search(possibles.begin(), possibles.end(), edgelist[i].second)) {
                            if (cliqueMembers.find(edgelist[i].first) == cliqueMembers.end())
                                continue;
                            neighbor1 = gr.neighbors(edgelist[i].second);
                            cliqueMembers.insert(edgelist[i].second);
                        } else {
                            disregard[i] = true;
                            continue;
                        }
                        disregard[i] = true;
                        possiblesNext.clear();
                        std::set_intersection(possibles.begin(), possibles.end(), neighbor1.begin(), neighbor1.end(), std::back_inserter(possiblesNext));
                        
                        // Swap possibles storage around
                        possibles.swap(possiblesNext);
                        
                        break;
                    }
                }
            }
            members = std::vector<int>(cliqueMembers.begin(), cliqueMembers.end());
        }
        for (itvec = members.begin(); itvec != members.end(); ++itvec) {
            cliques[(*itvec) - 1] = cliqueID;
        }
        cliqueID++;
        cliqueMembers.clear();
        gr.deleteNodes(members);
        nEdges = gr.nEdges();
    }
    
    if (gr.nNodes() != 0) {
        for (i = 0; i < cliques.size(); ++i) {
            if (cliques[i] == 0) {
                cliques[i] = cliqueID;
                cliqueID++;
            }
        }
    }
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

//[[Rcpp::export]]
DataFrame groupNeighbors(IntegerVector down, IntegerVector up, IntegerVector groups, IntegerVector order) {
    std::deque<int> groupOut, neighborOut;
    int currentGroup = groups[order[0]];
    int i, groupInd;
    std::set<int> neighborSet;
    
    for (i = 0; i < order.size(); ++i) {
        if (currentGroup != groups[order[i]]) {
            groupOut.insert(groupOut.end(), neighborSet.size(), currentGroup);
            neighborOut.insert(neighborOut.end(), neighborSet.begin(), neighborSet.end());
            neighborSet.clear();
            currentGroup = groups[order[i]];
        }
        groupInd = down[order[i]];
        if (groupInd != -1) {
            neighborSet.insert(groups[groupInd]);
        }
        groupInd = up[order[i]];
        if (groupInd != -1) {
            neighborSet.insert(groups[groupInd]);
        }
    }
    groupOut.insert(groupOut.end(), neighborSet.size(), currentGroup);
    neighborOut.insert(neighborOut.end(), neighborSet.begin(), neighborSet.end());
    
    return DataFrame::create(
        Named("group") = wrap(groupOut),
        Named("neighbor") = wrap(neighborOut)
    );
}

//[[Rcpp::export]]
DataFrame mergeGroupsByNeighbors(List GOI, DataFrame lookup) {
    int i, j, k, l;
    IntegerVector group;
    std::set<int> OGset;
    std::deque<int> pair1, pair2;
    bool pairFound;
    std::map< int, std::vector<int> > lookupMap;
    
    IntegerVector OG = lookup["OG"];
    IntegerVector NG = lookup["NG"];
    
    for (i = 0; i < OG.size(); ++i) {
        lookupMap[OG[i]].push_back(NG[i]);
    }
    
    // Start of horrible nesting - sorry
    for (i = 0; i < GOI.size(); ++i) {
        R_CheckUserInterrupt();
        
        group = GOI[i];
        
        for (j = 0; j < group.size() - 1; ++j) {
            OGset = std::set<int>(lookupMap[group[j]].begin(), lookupMap[group[j]].end());
            pairFound = false;
            
            for (k = j + 1; k < group.size(); ++k) {
                if (pairFound) break;
                
                for (l = 0; l < lookupMap[group[k]].size(); ++l) {
                    if (OGset.find(lookupMap[group[k]][l]) != OGset.end()) {
                        pair1.push_back(group[j]);
                        pair2.push_back(group[k]);
                        pairFound = true;
                        break;
                    }
                }
            }
        }
    }
    // end of horrible nesting...
    
    return DataFrame::create(
        Named("V1") = wrap(pair1),
        Named("V2") = wrap(pair2)
    );
}