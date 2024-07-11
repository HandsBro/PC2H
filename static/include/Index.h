#ifndef INDEX_H_
#define INDEX_H_
//--------------------------------------------------------------------------------------------------
#include "Graph.h"
#include "Timer.h"
//--------------------------------------------------------------------------------------------------
#include <map>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <string.h>
//--------------------------------------------------------------------------------------------------
#define NUM_THREADS 2 //for BFS only.
//--------------------------------------------------------------------------------------------------
struct queue_unit{
	int node;
	double score;
  int sum;
  queue_unit(int a, double sdg, int sm) : node(a), score(sdg), sum(sm) {}
  // queue_unit(int a, double sdg) : node(a), score(sdg), centrality(0) {}
  // queue_unit(int a, double sdg) : node(a), score(sdg) {}
	bool operator<(const queue_unit& a) const{
    // return score < a.score;//randomness makes the worst case slower.
    if (score == a.score) return sum < a.sum;
		else if(score > 0 && a.score>0) return score < a.score;
    else if(score > 0 && a.score==0) return false;
    else if(score == 0 && a.score>0) return true; 
    return false;
	}
};
// struct cmp{
//   bool operator()(queue_unit a,queue_unit b) { 
//   if (a.score == b.score) return a.centrality > b.centrality; 
//   return a.score < b.score; 
//   }
// };

// Hash function for a pair of vectors
struct VectorPairHash {
    std::size_t operator()(const std::pair<std::vector<unsigned>, std::vector<unsigned>>& pair) const {
        std::size_t seed = 0;
        for (const auto& v : pair.first) {
            seed ^= std::hash<int>{}(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        for (const auto& v : pair.second) {
            seed ^= std::hash<int>{}(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

// Equality function for a pair of vectors
struct VectorPairEqual {
    bool operator()(const std::pair<std::vector<unsigned>, std::vector<unsigned>>& lhs, const std::pair<std::vector<unsigned>, std::vector<unsigned>>& rhs) const {
        return lhs.first == rhs.first && lhs.second == rhs.second;
    }
};

class Index {
private:
  Graph *g;
  unsigned root;
  unsigned n_;
  unsigned Ctr;
  // int ranking;
  int labelsize;
  int labelsize_c;
  int max_level;
  int reachable;
  double debug_time = 0.0;
  double debug_time2 = 0.0;
  bool puredegree;
  int testsize;

  std::vector<unsigned> level;
  std::vector<unsigned> level2;
  std::vector<unsigned> level_c;
  std::vector<unsigned> level2_c;
  std::vector<std::vector<unsigned> > levelNodes; //classify nodes based on the levels.
  std::vector<double> inscores;
  std::vector<double> outscores;
  // std::vector<int> approximated_scores;
  std::vector<unsigned> calculated_order;
  std::priority_queue<queue_unit> pq_node_score;
  std::priority_queue<queue_unit> pq_node_score_leaves;
  std::priority_queue<queue_unit> pq_node_score_heavy;

  // std::priority_queue<queue_unit, std::vector<queue_unit>, cmp> pq_node_score;

  std::vector<std::vector<unsigned> > IN;
  std::vector<std::vector<unsigned> > OUT; 

  std::vector<std::vector<unsigned> > CIN;
  std::vector<std::vector<unsigned> > COUT;
  std::vector<unsigned> c_roots; 
  
  std::vector<int> visited;
  std::vector<int> redundant;
  std::vector<unsigned> * in_degree; 
  std::vector<unsigned> masked;
  std::vector<unsigned> maskedindex;

  //compression
  // std::unordered_map<std::pair<std::vector<unsigned>, std::vector<unsigned>>, std::vector<unsigned>, VectorPairHash, VectorPairEqual> vectorGroups;
  std::unordered_map<std::pair<std::vector<unsigned>, std::vector<unsigned>>, std::unordered_set<unsigned>, VectorPairHash, VectorPairEqual> vectorGroups;
  std::unordered_map<unsigned, std::unordered_set<unsigned>> ref_clusters;
  std::pair<std::vector<unsigned>, std::vector<unsigned>> vectors;//key
  std::vector<std::vector<unsigned> > cpd;
  std::vector<std::vector<unsigned> > cnb;
  std::vector<std::vector<unsigned> > cnb2;
  std::vector<unsigned> reduction_ref;
  std::vector<unsigned> origianl_ref;

  // std::vector<unsigned> querypointer;

  std::vector<unsigned> cin_degree;
  std::vector<unsigned> cout_degree;
  
  //delete
  std::vector<unsigned> delete_order;//insert in the reverse order.

  //query
  std::vector<unsigned> queriesbuildin;

  //ranks
  std::queue<unsigned> ranks;


public:
  /// constructor
  Index(Graph *g, unsigned rt);

  /// destructor
  ~Index();

  ///public
  unsigned PI();
  unsigned PI2();

  void out_score();
  void in_score();
  void calculate_order();
  void calculate_order_c();
  unsigned levelFilter(unsigned node);
  unsigned levelFilter_root(unsigned node);
  // unsigned levelFilter2(unsigned node);
  void _2hop();
  void _2hop_c();
  void output(std::string indexfile);
  void express();
  void storage();
  bool reachability(unsigned & u, unsigned & v);
  bool reachability_c(unsigned *u, unsigned *v);
  double get_debugtime();
  void clear_dummy();
  bool mergeCheck(unsigned *u, unsigned *v);
  bool mergeCheck_log(unsigned & u, unsigned & v);
  bool mergeCheck_c(unsigned *u, unsigned *v);

  void outputGraph();
  void outputQuery(std::vector<std::pair<unsigned, unsigned> > * Queries);

  void masknodes(int masksize);
  void nodes_delete();

  int compression();
  void compress_edges();

  void vector_delete(std::vector<unsigned> * vec, unsigned target);
  size_t hashvalue(std::pair<std::vector<unsigned>, std::vector<unsigned>> pair);

  void query(std::vector<std::pair<unsigned, unsigned> > *queries);
  void reachtest();

  void BFS(unsigned node);
  void BFStest();

  // accessors
  inline Graph* get_graph() const {
    return g;
  }

  //argmax/min
  template<class ForwardIterator>
  inline size_t argmin(ForwardIterator first, ForwardIterator last)
  {
      return std::distance(first, std::min_element(first, last));
  }

  template<class ForwardIterator>
  inline size_t argmax(ForwardIterator first, ForwardIterator last)
  {
      return std::distance(first, std::max_element(first, last));
  }
};
//--------------------------------------------------------------------------------------------------
#endif /* INDEX_H_ */
