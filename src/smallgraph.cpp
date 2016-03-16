#include "smallgraph.h"
#include <Rcpp.h>

using namespace Rcpp;

Graph::Graph(int nNodes, DataFrame& edges) {
    valid_edge = std::vector<bool>(edges.nrows(), true);
    next_possible_edge = 0;
    int i;
    for (i = 0; i < nNodes; ++i) {
        node_map.insert(std::make_pair(i + 1, std::vector<int>()));
    }
    from = edges["from"];
    to = edges["to"];
    edgecount = from.size();
    for (i = 0; i < from.size(); ++i) {
        node_map[from[i]].push_back(to[i]);
        node_map[to[i]].push_back(from[i]);
    }
    for (std::map< int, std::vector<int> >::iterator it = node_map.begin(); it != node_map.end(); ++it) {
        std::sort(it->second.begin(), it->second.end());
    }
}

void Graph::deleteNodes(std::vector<int> ids) {
    int i;
    std::map< int, std::vector<int> >::iterator itmap;
    std::vector<int>::iterator itvec;
    std::set<int>::iterator itset;
    std::sort(ids.begin(), ids.end());
    std::set<int> connectedNodes;
    IdInList id_in_list(ids);
    
    for (i = 0; i < ids.size(); ++i) {
        itmap = node_map.find(ids[i]);
        if (itmap != node_map.end()) {
            connectedNodes.insert(itmap->second.begin(), itmap->second.end());
            node_map.erase(itmap);
        }
    }
    for (i = 0; i < ids.size(); ++i) {
        itset = connectedNodes.find(ids[i]);
        if (itset != connectedNodes.end()) {
            connectedNodes.erase(itset);
        }
    }
    for (itset = connectedNodes.begin(); itset != connectedNodes.end(); ++itset) {
        itmap = node_map.find(*itset);
        if (itmap != node_map.end()) {
            itmap->second.erase(std::remove_if(itmap->second.begin(), itmap->second.end(), id_in_list), itmap->second.end());
        }
    }
    for (i = next_possible_edge; i < valid_edge.size(); ++i) {
        if (valid_edge[i]) {
            if (id_in_list(from[i]) || id_in_list(to[i])) {
                valid_edge[i] = false;
                --edgecount;
            }
        }
    }
}

std::vector<int> Graph::neighbors(int id) {
    std::map<int, std::vector<int> >::iterator itmap;
    itmap = node_map.find(id);
    if (itmap == node_map.end()) {
        stop("id not member of graph");
    }
    return itmap->second;
}

std::pair<int, int> Graph::firstEdge() {
    if (edgecount == 0) {
        stop("No edges in graph");
    }
    int i;
    for (i = next_possible_edge; i < valid_edge.size(); ++i) {
        if (valid_edge[i]) {
            next_possible_edge = i;
            break;
        }
    }
    return std::pair<int, int>(from[i], to[i]);
}

std::vector< std::pair<int, int> > Graph::completeTriangle(std::pair<int, int> edge) {
    std::set<int> commonNodes;
    std::vector< std::pair<int, int> > edges;
    std::map<int, std::vector<int> >::iterator itmap1, itmap2;
    std::set<int>::iterator itset;
    int i;
    
    itmap1 = node_map.find(edge.first);
    if (itmap1 == node_map.end()) {
        stop("Node not member of graph");
    }
    itmap2 = node_map.find(edge.second);
    if (itmap1 == node_map.end()) {
        stop("Node not member of graph");
    }
    std::set_intersection(itmap1->second.begin(), itmap1->second.end(), 
                          itmap2->second.begin(), itmap2->second.end(), 
                          std::inserter(commonNodes, commonNodes.begin()));
    commonNodes.insert(edge.first);
    commonNodes.insert(edge.second);
    for (i = next_possible_edge; i < valid_edge.size(); ++i) {
        if (valid_edge[i]) {
            itset = commonNodes.find(from[i]);
            if (itset == commonNodes.end()) continue;
            itset = commonNodes.find(to[i]);
            if (itset == commonNodes.end()) continue;
            edges.push_back(std::pair<int, int>(from[i], to[i]));
        }
    }
    return edges;
}

bool Graph::isComplete() {
    std::map< int, std::vector<int> >::iterator itmap;
    int maxNeighbors = node_map.size() - 1;
    for (itmap = node_map.begin(); itmap != node_map.end(); ++itmap) {
        if (itmap->second.size() != maxNeighbors) {
            return false;
        }
    }
    return true;
}

int Graph::nNodes() {
    return node_map.size();
}

std::vector<int> Graph::nodeIds() {
    std::vector<int> ids;
    for (std::map<int, std::vector<int> >::iterator it = node_map.begin(); it != node_map.end(); ++it) {
        ids.push_back(it->first);
    }
    return ids;
}

int Graph::nEdges() {
    return edgecount;
}

void Graph::show() {
    Rcout << "#Edges: " << edgecount << std::endl;
    for (int i = 0; i < valid_edge.size(); ++i) {
        if (valid_edge[i]) {
            Rcout << i << ": " << from[i] << ", " << to[i] << std::endl;
        }
    }
    Rcout << "#Nodes: " << node_map.size() << std::endl;
    std::map< int, std::vector<int> >::iterator itmap;
    for (itmap = node_map.begin(); itmap != node_map.end(); ++itmap) {
        Rcout << itmap->first << ": ";
        for (int k = 0; k < itmap->second.size(); ++k) {
            Rcout << itmap->second[k];
            if (k != itmap->second.size()-1) Rcout << ", ";
        }
        Rcout << std::endl;
    }
}