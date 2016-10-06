#ifndef SMALLGRAPH
#define SMALLGRAPH

#include <Rcpp.h>

using namespace Rcpp;

class IdInList {
    std::vector<int> * list;
    
public:
    IdInList (std::vector<int> & idlist) {
        list = &idlist;
    }
    
    bool operator()(int id) const {
        return std::binary_search(list->begin(), list->end(), id);
    }
};

class Graph{
public:
    Graph(int nNodes, DataFrame& edges);
    void deleteNodes(std::vector<int> ids);
    void deleteEdges(std::vector<int> ids);
    std::vector<int> neighbors(int id);
    int degree(int id);
    std::pair<int, int> firstEdge();
    std::vector< std::pair<int, int> > completeTriangle(std::pair<int, int> edge);
    bool isComplete();
    int nNodes();
    std::vector<int> nodeIds();
    std::vector<int> incidentEdges(int id);
    int nEdges();
    std::vector<int> edgeIds();
    std::pair<int, int> getEdge(int id);
    
    void show();
private:
    std::map<int, std::vector<int> > node_map;
    IntegerVector from;
    IntegerVector to;
    std::vector<bool> valid_edge;
    int next_possible_edge;
    int edgecount;
};

#endif
