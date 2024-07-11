#include "Graph.h"
#include <assert.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>
#include <string.h>
#include <algorithm>
//--------------------------------------------------------------------------------------------------
/**
 * Constructor
 * @param filename file containing adjacency lists
 */
//--------------------------------------------------------------------------------------------------
Graph::Graph(const std::string& filename) :
    n(0), m(0), both(0) {
  std::string line;
  std::ifstream stream(filename.c_str());
  if (!stream.eof()) {
    getline(stream, line);  //first line in the file provides information about #edge and #node
    std::istringstream iss(line, std::istringstream::in);
    iss >> n >> m; //m is useless in thie file, and only n will be accessed and used.
    nb = std::vector<std::vector<unsigned> >(n);  //adjacent list
    pd = std::vector<std::vector<unsigned> >(n);  //in list
    deg = std::vector<unsigned>(n, 0);
    // for (unsigned i = 0; i < n; ++i)
    //   deg[i] = 0;
    indeg = std::vector<unsigned>(n, 0);
    // leaves = std::vector<unsigned>(n, 0);
    unsigned s, t, d;
    std::vector<bool> is_root(n, true);
    getline(stream, line);
    while (stream.good()) {
      iss.clear();
      iss.str(line);
      iss >> s;
      d = 0;
      while (iss >> t) {  //file goes in this style: s t1 t2 t3 ...
        // std::cout << s << " -> " << t << std::endl;//no output
        nb[s].push_back(t);
        pd[t].push_back(s);
        ++d;
        ++indeg[t];
        is_root[t] = false;
      }
      // std::cout << (s < n) << n << std::endl;
      assert(s < n);
      deg[s] = d;
      getline(stream, line);
    }
    for (unsigned i = 0; i < n; ++i) {
      if (!deg[i]) {
        leaves.push_back(i);
      }
      if (is_root[i]) {
        roots.push_back(i);
      }
      if(!deg[i] && is_root[i]) {
        // std::cout<< "Isolated node ...\n";
        both++;
      }
    }
  }
}
//--------------------------------------------------------------------------------------------------
/**
 * Destructor
 */
//--------------------------------------------------------------------------------------------------
Graph::~Graph() {
}
//--------------------------------------------------------------------------------------------------
/**
 * Accessor for adjacency lists
 *
 * @param v vertex
 * @return adjacency list (std::vector)
 */
//--------------------------------------------------------------------------------------------------
// const std::vector<unsigned>* Graph::get_neighbors(unsigned v) const {
std::vector<unsigned>* Graph::get_neighbors(unsigned v) {
  return &nb[v];
}
//--------------------------------------------------------------------------------------------------
/**
 * Accessor for root vertices
 *
 * @return roots of the graph (vertices with zero indegree)
 */
//--------------------------------------------------------------------------------------------------
std::vector<unsigned>* Graph::get_roots() {
  return &roots;
}

int Graph::get_both() {
  return both;
}

//--------------------------------------------------------------------------------------------------
void Graph::clear_nb(unsigned v){
  nb[v].clear();
}

void Graph::clear_pd(unsigned v){
  pd[v].clear();
}
//--------------------------------------------------------------------------------------------------
void Graph::remove_nb(unsigned v, unsigned target){
  // std::cout<< nb[v].size() << std::endl;
  std::vector<unsigned>::iterator position = std::find(nb[v].begin(), nb[v].end(), target);
  if ( position!= nb[v].end()) nb[v].erase(position);
  // std::cout<< "after: " << nb[v].size() << std::endl;
}

void Graph::remove_pd(unsigned v, unsigned target){
  // std::cout<< nb[v].size() << std::endl;
  std::vector<unsigned>::iterator position = std::find(pd[v].begin(), pd[v].end(), target);
  if ( position!= pd[v].end()) pd[v].erase(position);
  // std::cout<< "after: " << nb[v].size() << std::endl;
}
//--------------------------------------------------------------------------------------------------
void Graph::add_root(unsigned node){
  roots.push_back(node);
}
