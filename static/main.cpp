#include "Graph.h"
#include "Timer.h"
#include "Index.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <vector>
#include <set>
#include <unordered_set>
#include <time.h>
#include <chrono>
// #include <unistd.h>

//read in queries
void read_queries(const std::string& query_file,
    std::vector<std::pair<unsigned, unsigned> > *queries) {
  std::ifstream qf(query_file.c_str(), std::ios::in);
  if (qf.is_open()) {
    std::string line;
    unsigned s, t;
    if (!qf.eof()) {
      while (qf.good()) {
        getline(qf, line);
        if (line.length()) {
          std::istringstream iss(line, std::istringstream::in);
          iss >> s;
          iss >> t;
          queries->push_back(std::make_pair(s, t));
        }
      }
    }
  }
}

int main(int argc, char *argv[]) {
  std::string graph_file = "", query_file = "", IndexFileStem = "";
  for (int i = 0; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-g") {
      graph_file = argv[++i];
    }
    else if (arg == "-q") {
      query_file = argv[++i];
    }
  }

  std::cout << "graph file: " << graph_file << std::endl;
  std::cout << "query file: " << query_file << std::endl;
  Graph *g_ = new Graph(graph_file);
  std::vector<std::pair<unsigned, unsigned> > queries;
  read_queries(query_file, &queries);

  unsigned rt = g_->num_nodes()-1;
  Index id(g_,rt);
  id.PI2();
  id.PI();
  std::cout << "Init Over ... " << std::endl;
  std::cout<<"Compressing...\n";
  Timer CP; CP.start();
  id.compression();
  id.compress_edges();
  double compress_time = CP.stop();
  id.in_score();
  id.out_score();
  std::cout<<"Indexing...\n";
  // Timer t1; t1.start();
  auto start = std::chrono::high_resolution_clock::now();
  id.calculate_order();
  id._2hop();
  auto stop = std::chrono::high_resolution_clock::now();
  double constructiontime = std::chrono::duration<double,std::milli>(stop-start).count();
  id.query(&queries);
  id.reachtest();
  std::cout << "Construction: (" << constructiontime << " ms)" << std::endl;
  std::cout << "Compressed time: (" << compress_time << " ms)" << std::endl;
  id.storage();
  std::cout<<"Making...\n";
  id.masknodes(10000);
  Timer t2; t2.start();
  id.nodes_delete();
  double Maintain = t2.stop();
  std::cout << "Maintain time: (" << Maintain << " ms)" << std::endl;

  delete g_;
  std::cout << "Done.\n";
  return 0;
  
}