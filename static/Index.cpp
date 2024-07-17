//! PC2H
#include "Index.h"
#include <assert.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string.h>
#include <queue>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <chrono>
#include <numeric>
#include <random>
//--------------------------------------------------------------------------------------------------
//! Constructor.
/*!
 * \param g pointer to graph
 * \param root dummy root of the graph.
 */
Index::Index(Graph* g_, unsigned rt){
  g = g_;
  n_ = g->num_nodes();
  root = rt;//not used
  labelsize = 0;
  labelsize_c = 0;
  reachable=0;
  levelNodes = std::vector<std::vector<unsigned> >(n_);
  puredegree = false;

  //level approximation.
  level = std::vector<unsigned>(n_,0);
  level2 = std::vector<unsigned>(n_,0);
  
  inscores = std::vector<double>(n_,0.00);
  outscores = std::vector<double>(n_,0.00);
  // approximated_scores = std::vector<int>(n_,0);
  // calculated_order 

  IN = std::vector<std::vector<unsigned> >(n_);
  OUT = std::vector<std::vector<unsigned> >(n_);

  visited = std::vector<int>(n_,-2);
  redundant = std::vector<int>(n_,0);
  in_degree = g->get_indegrees();

  reduction_ref = std::vector<unsigned>(n_);
  origianl_ref = std::vector<unsigned>(n_);
}
//--------------------------------------------------------------------------------------------------
//! Destructor.
Index::~Index() {
}
//--------------------------------------------------------------------------------------------------
unsigned Index::levelFilter(unsigned node){ 
  if (level[node] != 0){
      return 0;
  }
  // const std::vector<unsigned> *nb = g->get_neighbors(node);  //get the nb list of node.
  std::vector<unsigned> *nb = g->get_predecessors(node);  //get the nb list of node.
  for (std::vector<unsigned>::const_iterator it = nb->begin(); it != nb->end(); ++it) {
    levelFilter(*it);
  }

  if (nb->size() == 0){
      level[node] = 1;
      levelNodes[1].push_back(node);
      return 0;
  }
  else{
    unsigned maxChildren = 0;

    for (unsigned i=0; i<nb->size(); i++){
        if (level[(*nb)[i]] > maxChildren) {
            maxChildren = level[(*nb)[i]];
        }
    }
    level[node] = maxChildren + 1;
    levelNodes[maxChildren + 1].push_back(node);
    return 0;
  }
}

unsigned Index::PI2(){//general case of level filters
  // std::vector<unsigned> * indegree = g->get_indegrees();//for degree based
  // std::vector<unsigned> * degree = g->get_degrees();
  unsigned useless = 0;
  // const std::vector<unsigned> *nb = g->get_roots();  //get the nb list of node:ignore dummuy root.
  const std::vector<unsigned> *nb = g->get_leaves();  //get the nb list of node:ignore dummuy root.
  // Timer t1; t1.start();
  for (unsigned i=0; i<nb->size(); ++i) {
    // std::cout<< i << std::endl;
    useless = levelFilter((*nb)[i]);
  }
  max_level = *max_element(level.begin(), level.end());
  if (max_level==2) puredegree=true;

  std::cout<< "Max Level: " << max_level << std::endl;
  if((int)(n_- g->get_roots()->size() - nb->size() + g->get_both())<1000000) puredegree=true;
  return useless;
}

unsigned Index::levelFilter_root(unsigned node){
  // std::cout<<"Reaching node: "<<node<< "levels: " << level[node] <<std::endl; //debug process
  if (level2[node] != 0){
      return 0;
  }
  const std::vector<unsigned> *nb = g->get_neighbors(node);  //get the nb list of node.
  // std::vector<unsigned> *nb = g->get_predecessors(node);  
  for (std::vector<unsigned>::const_iterator it = nb->begin(); it != nb->end(); ++it) {
    levelFilter_root(*it);
  }

  if (nb->size() == 0){
      level2[node] = 1;
      // levelNodes[1].push_back(node);
      return 0;
  }
  else{
    unsigned maxChildren = 0;

    for (unsigned i=0; i<nb->size(); i++){
        if (level2[(*nb)[i]] > maxChildren) {
            maxChildren = level2[(*nb)[i]];
        }
    }
    level2[node] = maxChildren + 1;
    // levelNodes[maxChildren + 1].push_back(node);
    return 0;
  }
}

unsigned Index::PI(){//from roots
  // const std::vector<unsigned> *roots = g->get_roots(); 
  // visited = std::vector<int>(n_,0);
  unsigned useless = 0;
  const std::vector<unsigned> *nb = g->get_roots();  //get the nb list of node:ignore dummuy root.
  for (unsigned i=0; i<nb->size(); ++i) {
    // std::cout<< i << std::endl;
    useless = levelFilter_root((*nb)[i]);
  }; 
  // std:: cout<< level2[1] << "@@@@@@@@@\n";
  return useless;
}

void Index::out_score(){
  // std::vector<unsigned>* deg = g->get_indegrees();
  // double limit = 1e20;
  double out_sum = 0;
  for (int i = max_level; i>=1; i--){
    std::vector<unsigned> levelset = levelNodes[i];
    if (i==max_level){//set to 0
      for (std::vector<unsigned>::iterator it=levelset.begin(); it!=levelset.end(); ++it){
        outscores[*it] = 0.00;
      }
    }

    else{
      for (std::vector<unsigned>::iterator it=levelset.begin(); it!=levelset.end(); ++it){
        const std::vector<unsigned> *nb = g->get_neighbors(*it); 
        for (std::vector<unsigned>::const_iterator it2=nb->begin(); it2!=nb->end();++it2){
          outscores[*it] += outscores[*it2]+1;
          // if(outscores[*it]>=limit) outscores[*it]=limit;
          // if(level[*it2]-level[*it]>1) jummping_edges++;
          // total_edges++;
          // outscores[*it] = outscores[*it] + (outscores[*it2] + 1)/(*deg)[*it2];
        }
        out_sum += outscores[*it];
      }
    }
  }
  // std::cout<< "total edges: " << total_edges << "; jumping edges: " << jummping_edges << "; ratio: " << jummping_edges/total_edges << std::endl;
  // std::cout<< "total nodes (without roots): " << n_- g->get_leaves()->size()<< " ; out sum: " << out_sum << std::endl;
}

void Index::in_score(){
  // std::vector<unsigned> * indegree = g->get_indegrees();
  // std::vector<unsigned> * degree = g->get_degrees();
  // std::vector<unsigned>* deg = g->get_degrees();
  // double limit = 1e20;
  double in_sum = 0;
  for (int i = 1; i<=max_level; ++i){
    std::vector<unsigned> levelset = levelNodes[i];
    if (i==1){//set to 0
      for (std::vector<unsigned>::iterator it=levelset.begin(); it!=levelset.end(); ++it){
        inscores[*it] = 0.00;
      }
    }
    else{
      for (std::vector<unsigned>::iterator it=levelset.begin(); it!=levelset.end(); ++it){
        const std::vector<unsigned> *pd = g->get_predecessors(*it); 
        for (std::vector<unsigned>::const_iterator it2=pd->begin(); it2!=pd->end();++it2){
          inscores[*it] += inscores[*it2]+1;
          // if(inscores[*it]>=limit) inscores[*it]=limit;
          // inscores[*it] = inscores[*it] + (inscores[*it2] + 1)/(*deg)[*it2];
        }
        in_sum += inscores[*it];
      }
    }
  }
  // std::cout<< "In socre[1]: " << inscores[1] << "; indeg: " << (*indegree)[1] << "; outdeg: " << (*degree)[1] << std::endl;
  // std::cout<< "total nodes (without roots): " << n_- g->get_leaves()->size()<< " ; in sum: " << in_sum << std::endl;
}

void Index::clear_dummy(){
  const std::vector<unsigned> *nb = g->get_neighbors(root); 
  visited[root] = -1;
  for (std::vector<unsigned>::const_iterator it=nb->begin(); it!=nb->end();++it){
    // if (outscores[*it] > 0) inscores[*it] = 0.0;
    (*in_degree)[*it]=0;
  }
}

void Index::BFS(unsigned node){
  std::cout << node << " @.@ " << std::endl;
  // std::vector<int> color(n_,0);
  std::queue<unsigned> q;
  int meets=0;
  std::vector<unsigned> meetlist;
  std::set<unsigned> meetset;
  // q.push(node);
  // while(!q.empty()){
  //   int temp = q.front();
  //   q.pop();
  //   if(visited[temp]==(int)node || visited[temp]==-1) continue;
  //   visited[temp]=node;
  //   // IN[temp].push_back(1);
  //   const std::vector<unsigned> *nb = g->get_neighbors(temp);
  //   for(auto it: *nb) q.push(it);
  // }
  visited[node]=-2;
  q.push(node);
  // double getp_time = 0.0;
  while(!q.empty()){
    int temp = (q.front());
    q.pop();
    if(visited[temp]==(int)node) continue;
    meets++;
    // meetlist.push_back(temp);
    // meetset.insert(reduction_ref[temp]);
    visited[temp]=node;
    
    std::vector<unsigned> *pd = g->get_predecessors(temp);
    // Timer t1;
    // t1.start();
    
    for(auto & it: (*pd)) q.push(it);
    // double bfstime = t1.stop();
    // getp_time += bfstime;
  }
  // visited[node]=-1;
  // std::cout << "push time is " << getp_time << std::endl;
  // std::cout << "BFS test time is: " << bfstime << std::endl;
  std::cout << "meets is: " << meets << std::endl;
  // sort(meetlist.begin(), meetlist.end());
  // for(auto it: meetset){
  //   std::cout << it << std::endl;
  // }
}

void Index::BFStest(){
  // const std::vector<unsigned> *rt = g->get_roots(); 
  // std::cout << "roots size is: " << rt->size() << std::endl;
  Timer t1; t1.start();
  // for(auto it: *rt) BFS(it);
  // while(!ranks.empty()){
  BFS(717112);//G^C: 1897335 G: 1824817    shuffle test 2437517
  //citeseerx: 717112
  // ranks.pop();
  // break;
  // }
  double bfstime = t1.stop();
  std::cout << "BFS test time is: " << bfstime << std::endl;

}

void Index::calculate_order(){
  std::vector<unsigned> * indegree = g->get_indegrees();//for degree based
  std::vector<unsigned> * degree = g->get_degrees();
  // visited = std::vector<int>(n_,0);
  // int in2de, out2de;
  std::vector<unsigned> temp(vectorGroups.size(),0);
  if(puredegree){
    for (unsigned i=0; i<n_; ++i){
      // if(ref_clusters[reduction_ref[i]].size()>1) continue;
      // if(temp[reduction_ref[i]]==1) {
      //   visited[i]=-1;
      //   continue;
      // }
      // temp[reduction_ref[i]]=1;

      //! 1-hop version.
      double temp_score;
      // if(i==1) temp_score = 1e20;
      temp_score = (double)((*indegree)[i]) * (double)((*degree)[i]) / ((double)((*indegree)[i]) + (double)((*degree)[i]) + 1);

      queue_unit temp(i, temp_score,(*degree)[i]+(*indegree)[i]);//1-hop version
      // queue_unit temp(i, temp_score,in2de+out2de);
      pq_node_score.push(temp);
    }
  }
  else{
    for (unsigned i=0; i<n_; ++i){
      // if(temp[reduction_ref[i]]==1) {
      //   visited[i]=-1;
      //   continue;
      // }
      // temp[reduction_ref[i]]=1;

      double temp_score = (inscores[i]*outscores[i])/(inscores[i]+outscores[i]+1);
      // double temp_score = (inscores[i]*outscores[i]);
      // if (temp_score < 0) std::cout<< "OVERFLOW" << std::endl;

      queue_unit temp(i,temp_score,inscores[i]+outscores[i]);
      pq_node_score.push(temp);

    }
    
  }

  // std::ifstream file("./ranks/rank_citp_c.txt");
  // std::string line;
  // if (file.is_open()) {
  //     while (std::getline(file, line)) {
  //         if (!line.empty()) {
  //             // ranks.push(std::stoi(line));//1897335
  //             // ranks.push(origianl_ref[std::stoi(line)]); //1824817
  //             // ranks.push(1824817);
  //         }
  //     }
  //     file.close();
  // } else {
  //     std::cerr << "Unable to open Rank file" << std::endl;
  //     return;
  // }
  // std::cout << "Rank Reading Done; Build-in Size: " << pq_node_score.size() << std::endl;
}

void Index::_2hop(){
  //original graph
  std::vector<unsigned> * indegree = g->get_indegrees();
  std::vector<unsigned> * degree = g->get_degrees();
  // std::vector<int> output_visited(n_,0);
  // visited = std::vector<int>(n_,0);
  // std::vector<unsigned> q;
  std::queue<unsigned > q;
  int i = 0;
  unsigned currentNode;
  
  // double uptime=0.00;
  // double downtime=0.00;
  // double pq_time=0.00;
  // int up_check_time=0, down_check_time=0;
  // Timer t1;
  // int up_success_check_time=0, down_success_check_time=0;
  double temp_score;
  int temp_sum;

  // std::ofstream ranking("./ranks/rank_citp_Gtest.txt");
  while (!pq_node_score.empty()){
  // while (!ranks.empty()){
    // t1.start();/
    while(1){
      int top_node = pq_node_score.top().node;
      double top_score = pq_node_score.top().score;
      int top_sum = pq_node_score.top().sum;

      if(!puredegree) 
      {
        temp_score = (inscores[top_node]*outscores[top_node])/(inscores[top_node]+outscores[top_node]+1);
        temp_sum = inscores[top_node]+outscores[top_node];
      }
      // if(!puredegree) temp_score = (inscores[top_node]*outscores[top_node]);
      else {
        // if(top_node==1) temp_score = 1e20;
        temp_score = (double)((*indegree)[top_node]) * (double)((*degree)[top_node]) / (double)((*indegree)[top_node] + (*degree)[top_node] + 1) ;//new
        temp_sum = (*indegree)[top_node] + (*degree)[top_node];
      }

      if (temp_score == top_score && temp_sum == top_sum){
        currentNode = top_node;
        pq_node_score.pop();
        break;
      }
      pq_node_score.pop();

      if(!puredegree){
        queue_unit temp_struct(top_node, temp_score, inscores[top_node]+outscores[top_node]);
        pq_node_score.push(temp_struct);
      } 
      else{
        queue_unit temp_struct(top_node, temp_score, temp_sum);
        pq_node_score.push(temp_struct);
      } 
    }
    // pq_time+=t1.stop();

    // if(output_visited[reduction_ref[currentNode]]==0){
    //   output_visited[reduction_ref[currentNode]]=1;
    //   // ranking << reduction_ref[currentNode] << std::endl;
    //   ranking << currentNode << std::endl;//for c G.
    // }

    // currentNode = ranks.front();
    // ranks.pop();

    ++i;
    if(visited[currentNode]==-1) continue;
    if ((*indegree)[currentNode] + (*degree)[currentNode] == 0) {
      // IN[currentNode].emplace_back(i);
      // OUT[currentNode].emplace_back(i);
      break;
    } 
    // std::cout<< i << std::endl;

    q.push(currentNode);
    while(!q.empty()){
      unsigned temp = (q.front());
      q.pop();
      // down_check_time++;

      if((temp != currentNode)){
        inscores[temp] -= (inscores[currentNode]+1);
        if(inscores[temp]<0) std::cout<< "NEGATIVE inscore:" << inscores[temp] <<std::endl;
      }

      if (visited[temp]==i || visited[temp]==-1) {
        continue;//visited before
      }
      visited[temp]=i;

      // desc++;

      if (temp==currentNode) {
        // IN[reduction_ref[temp]].emplace_back(i);
        IN[temp].emplace_back(i);
        std::vector<unsigned> *nb = g->get_neighbors(temp);
        // unsigned* ptr = nb->data();
        // for (unsigned k=0; k<nb->size(); ++k){
        // for (unsigned k=0; k<cnb[temp].size(); ++k){
        // for(auto it=nb->begin(); it!=nb->end(); it++){
        for(auto & it: *nb){
          (*indegree)[it] -= 1;
          q.push(it);
          // q.push(cnb[temp][k]);
        }

      }
      else {
        // down_check_time++;
        // if (!mergeCheck_log(reduction_ref[currentNode], reduction_ref[temp])) {
        // t1.start();
        if (!mergeCheck_log(currentNode,temp)) {
        // if(true){
        // if (!mergeCheck(&currentNode,&temp)) {
          // IN[reduction_ref[temp]].emplace_back(i);
          IN[temp].emplace_back(i);
          labelsize += 1;

          // if(!puredegree){
          //   inscores[temp] -= (inscores[currentNode]+1);
          //   // if(inscores[temp]<0) std::cout<< "NEGATIVE inscore:" << inscores[temp] <<std::endl;
          // }

          std::vector<unsigned> *nb = g->get_neighbors(temp);
          // unsigned* ptr = nb->data();
          // for (unsigned k=0; k<nb->size(); ++k){
          // for (unsigned k=0; k<cnb[temp].size(); ++k){
          // for(auto it=nb->begin(); it!=nb->end(); it++){
          for(auto & it: *nb){
            q.push(it);
            // q.push(cnb[temp][k]);
          }
        }
        // downtime += t1.stop();
      }
    }
    //====
    visited[currentNode]=-2;
    q.push(currentNode);
    while(!q.empty()){
      // unsigned temp = q[z];
      unsigned temp = (q.front());
      q.pop();
      // up_check_time++;

      if(temp != currentNode){
        outscores[temp] -= (outscores[currentNode]+1);
        if(outscores[temp]<0) std::cout<< "NEGATIVE outscores:" << outscores[temp] <<std::endl;
      }

      if (visited[temp]==i || visited[temp]==-1) {
        continue;//visited before
      }

      visited[temp]=i;

      if (temp == currentNode) {
        // OUT[reduction_ref[temp]].emplace_back(i);
        OUT[temp].emplace_back(i);

        std::vector<unsigned> *pd = g->get_predecessors(temp);
        // unsigned* ptr = pd->data();
        // for (unsigned k=0; k<pd->size(); ++k){
        // for(auto it=pd->begin(); it!=pd->end(); it++){
        // for (unsigned k=0; k<cpd[temp].size(); ++k){
        for(auto & it: *pd){
          (*degree)[it] -= 1;//update degree.
          q.push(it);
          // q.push(cpd[temp][k]);
        }
      }

      else {
        // up_check_time++;
        // if (!mergeCheck_log(reduction_ref[temp], reduction_ref[currentNode])) {
        // t1.start();
        if (!mergeCheck_log(temp, currentNode)) {
        // if(true){
        // if (!mergeCheck(&temp,&currentNode)) {
          // up_success_check_time++;
          // OUT[reduction_ref[temp]].emplace_back(i);
          OUT[temp].emplace_back(i);
          labelsize +=1;
          
          std::vector<unsigned> *pd = g->get_predecessors(temp);
          // unsigned* ptr = pd->data();
          // for (unsigned k=0; k<pd->size(); ++k){
          // for(auto it=pd->begin(); it!=pd->end(); it++){
          // for (unsigned k=0; k<cpd[temp].size(); ++k){
          for(auto & it: *pd){
            // q.push((*pd)[k]);
            q.push(it);
            // q.push(cpd[temp][k]);
          }
        }
        // uptime += t1.stop();
      }
    }

    visited[currentNode] = -1;//set to -1 to permanently mask.

  }
  // std::cout<< "Forward BFS: " << downtime << "; Back BFS: " << uptime << std::endl;
  // printf("Down total: %d; Up total: %d.\n",down_check_time,up_check_time);
  // printf("Down total: %d; Positive total: %d; Up total: %d; positive total: %d.\n",down_check_time, down_success_check_time,up_check_time, up_success_check_time);
  // printf("pq time is: %f\n", pq_time);

  // ranking.close();
}

void Index::output(std::string st){
  std::string s1 = st + ".label";
  std::string s2 = st + ".index";
  FILE *file = fopen(s1.c_str(),"w");//label
  FILE *ifile = fopen(s2.c_str(),"w");//index


  int x = 0;
  for ( unsigned i = 0 ; i < n_ ; ++i )
  {
    fprintf(ifile,"%d ",x);

    fprintf(file,"%d %d ",level[i], level2[i]); 
    x+=2;
    
    for ( unsigned j = 0 ; j <OUT[i].size() ; ++j )
    {
      fprintf(file,"%d ",OUT[i][j]);
      x++;
    }

    fprintf(file,"%d ",-1);
    x++;

    fprintf(ifile,"%d\n",x);

    fprintf(file,"%d %d ",level[i], level2[i]); 
    x+=2;
    
    for ( unsigned j = 0 ; j < IN[i].size() ; ++j )
    {
      fprintf(file,"%d ",IN[i][j]);
      x++;
    }

    fprintf(file,"%d\n",-2);
    x++;
  }

  fclose(file);
  fclose(ifile);
}

void Index::query(std::vector<std::pair<unsigned, unsigned> > *queries){
  testsize=1'000'000;//100000
  queriesbuildin.resize(2*testsize);
  for(int i=0; i<testsize; ++i){
    // queriesbuildin[i] = origianl_ref[reduction_ref[(*queries)[i].first]];
    // queriesbuildin[i+testsize] = origianl_ref[reduction_ref[(*queries)[i].second]]; 
    queriesbuildin[i] = (*queries)[i].first;
    queriesbuildin[i+testsize] = (*queries)[i].second; 
  }
}

void Index::reachtest(){
  reachable=0;
  auto start = std::chrono::high_resolution_clock::now();
  for(int k=0; k<testsize; ++k){
    unsigned u = queriesbuildin[k], v = queriesbuildin[k+testsize];
    if (level[u] >= level[v] || level2[u] <= level2[v]) {
      // std::cout<< 0 << std::endl;
      continue;
    }
    if (IN[v].size()==0 || OUT[u].size()==0) {
      // std::cout<< 0 << std::endl;
      continue; 
    }
    //merge check
    unsigned i=0, j=0;
    while(1){
      // std::cout<< i << ";;" << j << std::endl;
      if (OUT[u][i] == IN[v][j]) {
        reachable++;

        // std::cout<< 1 << std::endl;
        break;
      }
      if (OUT[u][i] > IN[v][j]){
        ++j;
        if (j==IN[v].size()) {
          // std::cout<< 0 << std::endl;
          break;
        }
      }
      if (OUT[u][i] < IN[v][j]){
        ++i;
        if (i==OUT[u].size()) {
          // std::cout<< 0 << std::endl;
          break;
        }
      }
    }    

    // bs
    // if (IN[v].size() <= OUT[u].size()){
    //   for(unsigned & it: IN[v]){
    //   // for(std::vector<unsigned>::iterator it = IN[v].begin(); it != IN[v].end(); ++it){
    //     if(it==*(std::lower_bound(OUT[u].begin(),OUT[u].end(), it))) {
    //       reachable++;
    //       break;
    //     }
    //   }
    // }
    // else{
    //   for(unsigned & it: OUT[u]){
    //   // for(std::vector<unsigned>::iterator it = OUT[u].begin(); it != OUT[u].end(); ++it){
    //     if(it==*(std::lower_bound(IN[v].begin(),IN[v].end(), it))){
    //       reachable++;
    //       break;
    //     } 
    //   }
    // }
  }
  auto stop = std::chrono::high_resolution_clock::now();
	double reachabilityTime = std::chrono::duration<double,std::milli>(stop-start).count();
  std::cout << "Query Time: (" << reachabilityTime << " ms)" <<std::endl;
}

bool Index::reachability(unsigned &u, unsigned &v){
  // std::cout<< OUT[u].size()<< "; " << IN[v].size() << std::endl; //371764
  // if (u==v) return true;
  // std::cout<< level[*u]<< "; " << level[*v] << std::endl;
  if (level[u] >= level[v] || level2[u] <= level2[v]) return false;
  // if (l) return false;
  // return mergeCheck_log(reduction_ref[u],reduction_ref[v]);
  return mergeCheck_log(u,v);
  // return mergeCheck(&u,&v);
}

bool Index::reachability_c(unsigned *u, unsigned *v){
  if (*u==*v) {
    return true;  }
  if (level[*u] >= level[*v]) return false;
  if (level2[*u] <= level2[*v]) return false;
  // std::cout<< reduction_ref[*u] << "?->" << reduction_ref[*v] << std::endl;
  return mergeCheck_c(&reduction_ref[*u],&reduction_ref[*v]);
  // return mergeCheck_c(u,v);
}

bool Index::mergeCheck(unsigned *u, unsigned *v){
  debug_time2++;
  unsigned i=0, j=0;
  if (IN[*v].size()==0 || OUT[*u].size()==0) return false;
  // std::cout<< "u: " << *u << ", u size: " <<  OUT[*u].size() << "; v: " << *v << ", v size: " << IN[*v].size() << std::endl;
  // std::cout<< "level u:" << level2[*u] << ", u back: " <<  OUT[*u].back() << "; level[v]: " << level2[*v] << ", v back: " << IN[*v].back() << std::endl;
  debug_time+=2;
  while(1){
    // std::cout<< i << ";;" << j << std::endl;
    if (OUT[*u][i] == IN[*v][j]) {
      reachable++;
      // std::cout<< "u: " << *u << ", u size: " <<  OUT[*u].size() << "; v: " << *v << ", v size: " << IN[*v].size() << std::endl;
      // std::cout<< "labels accessed: " << local << std::endl;
      return true;
    }
    if (OUT[*u][i] > IN[*v][j]){
      ++j;
      debug_time++;
      if (j==IN[*v].size()) {
        // std::cout<< "labels accessed: " << local << std::endl;
        return false;
      }
    }
    if (OUT[*u][i] < IN[*v][j]){
      ++i;
      debug_time++;
      if (i==OUT[*u].size()) {
        // std::cout<< "labels accessed: " << local << std::endl;
        return false;
      }
    }
  }
}

bool Index::mergeCheck_log(unsigned &u, unsigned &v){
  // std::cout<< "u: " << *u << ", u size: " <<  OUT[*u].size() << "; v: " << *v << ", v size: " << IN[*v].size() << std::endl;
  if (IN[v].size()==0 || OUT[u].size()==0) return false;
  if (IN[v].size() <= OUT[u].size()){
    for(unsigned it: IN[v]){
    // for(std::vector<unsigned>::iterator it = IN[v].begin(); it != IN[v].end(); ++it){
      if(it==*(std::lower_bound(OUT[u].begin(),OUT[u].end(), it))) {
        // reachable++;
        return true;
      }
    }
    return false;
  }
  else{
    for(unsigned it: OUT[u]){
    // for(std::vector<unsigned>::iterator it = OUT[u].begin(); it != OUT[u].end(); ++it){
      if(it==*(std::lower_bound(IN[v].begin(),IN[v].end(), it))){
        // reachable++;
        return true;
      } 
    }
    return false;
  }
}

bool Index::mergeCheck_c(unsigned *u, unsigned *v){
  // std::cout<< OUT[*u].size() << "; " << IN[*v].size() << std::endl;
  unsigned i=0, j=0;
  if (CIN[*v].size()==0 || COUT[*u].size()==0) return false;
  while(1){
    if (COUT[*u][i] == CIN[*v][j]) {
      reachable++;
      return true;
    }
    if (COUT[*u][i] > CIN[*v][j]){
      ++j;
      if (j==CIN[*v].size()) return false;
    }
    if (COUT[*u][i] < CIN[*v][j]){
      ++i;
      if (i==COUT[*u].size()) return false;
    }
  }
}

void Index::storage(){
  float totalsize = labelsize * 4;
  // float totalsize_c = labelsize_c * 4;
  unsigned MAXOUT=0, MAXIN=0;//, MAXIN_id;
  for(unsigned i = 0; i<n_; ++i){
    if(IN[i].size()>MAXIN) {
      MAXIN = IN[i].size();
    }
    if(OUT[i].size()>MAXOUT) MAXOUT = OUT[i].size();

  }

  std::cout<< "Index size (byte): " << totalsize << std::endl;
  std::cout << "Total Reachable: " << reachable <<std::endl;
}

double Index::get_debugtime(){
  return debug_time;
}

void Index::express(){
  int c = 0;
  std::vector<unsigned> *nb;
  std::vector<unsigned> *pd;
  for (auto it = vectorGroups.begin(); it != vectorGroups.end(); ++it) {
      std::cout << "Vector 1: ";
      for (int val : it->first.first) {
          std::cout << val << " ";
      }
      std::cout << ", Vector 2: ";
      for (int val : it->first.second) {
          std::cout << val << " ";
      }
      std::cout << std::endl;
      std::cout << "Numbers: ";
      for (int num : it->second) {
          nb = g->get_neighbors(num);
          pd = g->get_predecessors(num);
          std::cout << num << ": level " << level[num] << "; nb&pd sizes: " << nb->size() << "," << pd->size() << "; ";
      }
      std::cout << std::endl;
      c++;

      // if (c==100) break;
  }
}

int Index::compression(){
  // int k=0;
  // for (unsigned i=0; i<n_; ++i){
  // std::unordered_set<unsigned>redundants;
  for (unsigned layer = max_level; layer>=1; layer--){
    for(unsigned i: levelNodes[layer]){
      // std::cout<< k++ << std::endl;
      
      //version 1: with reachability indices.
      std::vector<unsigned> * nb = g->get_neighbors(i);
      std::vector<unsigned> * pd = g->get_predecessors(i);
      // std::vector<unsigned> nb2(g->get_neighbors(i)->begin(),g->get_neighbors(i)->end());
      // std::vector<unsigned> pd2(g->get_predecessors(i)->begin(),g->get_predecessors(i)->end());

      // printf("Now at node: %d with nb&pd size: %lu,%lu\n", i, nb->size(), pd->size());
      //eliminate edges

      // for(auto t: nb2){//origianl by out(i)
      //   for(auto s: nb2){
      //     if(s==t) continue;
      //     if (mergeCheck(&s, &t)){//delete
      //       std::vector<unsigned>::iterator position = std::find(nb->begin(), nb->end(), t);
      //       if ( position!= nb->end()) {
      //         nb->erase(position);
      //         std::vector<unsigned> * pd_temp = g->get_predecessors(t);
      //         position = std::find(pd_temp->begin(), pd_temp->end(), i);
      //         pd_temp->erase(position);
      //         break;
      //       }
      //     }
      //   }
      // }

      // redundants.clear();
      // for(auto t: nb2){//check by min{o,i}
      //   std::vector<unsigned> * pdt = g->get_predecessors(t);
      //   // std::cout << l++ << "==> "<< "pd size: " << pdt->size() << "; nb size: " << nb->size() << std::endl;
      //   if(pdt->size()==1) continue;
      //   if(redundants.count(t)) continue;
      //   if(pdt->size()>=nb->size()){//nb smaller
      //     for(auto s: nb2){
      //       if(s==t) continue;
      //       if (mergeCheck(&s, &t)){//delete
      //         redundants.insert(t);
      //         std::vector<unsigned>::iterator position = std::find(nb->begin(), nb->end(), t);
      //         if ( position!= nb->end()) {
      //           nb->erase(position);
      //           std::vector<unsigned> * pd_temp = g->get_predecessors(t);
      //           position = std::find(pd_temp->begin(), pd_temp->end(), i);
      //           pd_temp->erase(position);
      //           break;
      //         }
      //       }
      //     }
      //   }
      //   else{//pd smaller
      //     for(auto s: *pdt){
      //       if(s==i) continue;
      //       if (mergeCheck(&i, &s)){//delete
      //         redundants.insert(t);
      //         std::vector<unsigned>::iterator position = std::find(nb->begin(), nb->end(), t);
      //         if ( position!= nb->end()) {
      //           nb->erase(position);
      //           std::vector<unsigned> * pd_temp = g->get_predecessors(t);
      //           position = std::find(pd_temp->begin(), pd_temp->end(), i);
      //           pd_temp->erase(position);
      //           break;
      //         }
      //       }      //     }
      //   }
      // }

      std::sort(nb->begin(), nb->end());// does not really influence the compression time.
      std::sort(pd->begin(), pd->end());
      vectors = std::make_pair(*nb, *pd);
      vectorGroups[vectors].insert(i);

      // vectorGroups[vectors].push_back(i);
    }
  }
  // printf("|V| before compressing: %u\n", n_);
  // printf("|V| after compressing: %lu\n", vectorGroups.size());


  // for (auto it = vectorGroups.begin(); it != vectorGroups.end(); ++it) {
  //   if((it->second).size()>1) std::cout<< "level: " << level[*(it->second.begin())]<<"; size: " << (it->second).size() << std::endl;
  // }

  return vectorGroups.size(); //for query construction
}

void Index::compress_edges(){
  // std::vector<unsigned> * indegree = g->get_indegrees();
  // std::vector<unsigned> * degree = g->get_degrees();
  unsigned c = 0;
  // unsigned cedges = 0;
  unsigned cedges2 = 0;
  int new_size = vectorGroups.size();
  // global
  // cpd = std::vector<std::vector<unsigned> >(new_size);//no need to copy when constructing.
  cnb = std::vector<std::vector<unsigned> >(new_size);
  // cnb2 = std::vector<std::vector<unsigned> >(n_);

  // cin_degree = std::vector<unsigned>(new_size);
  // cout_degree = std::vector<unsigned>(new_size);

  // CIN = std::vector<std::vector<unsigned> >(new_size);
  // COUT = std::vector<std::vector<unsigned> >(new_size);

  //! randomized id.
  // std::vector<unsigned> shufflebase(new_size);
  // std::iota(shufflebase.begin(), shufflebase.end(), 0);
  // std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
  // std::shuffle(shufflebase.begin(), shufflebase.end(), generator);
  
  // unsigned test_;
  //namespace: cluster order;
  for (auto it = vectorGroups.begin(); it != vectorGroups.end(); ++it) {
    // unsigned size_cnt=0; 
    // int newid = shufflebase[c];
    int newid = c;
      for (int num : it->second) {
        // size_cnt++;
        // if(size_cnt<it->second.size()) {
        //   // visited[num]=-1;
        //   redundant[num]=1;
        // }
        reduction_ref[num] = newid;
        origianl_ref[newid] = num;
        ref_clusters[newid].insert(num);
      }
      c++;
  }


  std::vector<unsigned> signs(c,-1);
  c=0;

  //! nb by cluster.
  for (auto it = vectorGroups.begin(); it != vectorGroups.end(); ++it) {
    // unsigned newid = shufflebase[c];
    unsigned newid = c;
    //nb
    for (auto val : it->first.first) {
      if(signs[reduction_ref[val]]!=newid){
        cnb[newid].push_back(reduction_ref[val]);
        signs[reduction_ref[val]]=newid;
      }
      
    }
    //pd
    // for (int val : it->first.second) {
    //   if(signs[reduction_ref[val]]!=newid){
    //     cpd[newid].push_back(reduction_ref[val]);
    //     signs[reduction_ref[val]]=newid;
    //   }
    // }

    // cedges += cpd[newid].size();
    cedges2 += cnb[newid].size();
    c++;//in the same order of addition.
  }

  //verification process: ensure each node in c has the same neighbour.
  // for (auto it = vectorGroups.begin(); it != vectorGroups.end(); ++it) {
  //   if(1){
  //     for (int num : it->second) {
  //       std::vector<unsigned> * nb = g->get_neighbors(num);
  //       std::vector<unsigned> * pd = g->get_predecessors(num);
  //       for(unsigned i=0; i<nb->size(); ++i){
  //         if((*nb)[i]!=it->first.first[i]) std::cout << "ERROR COMPRESS!!!!!!!\n";
  //       }

  //       for(unsigned i=0; i<pd->size(); ++i){
  //         if((*pd)[i]!=it->first.second[i]) std::cout << "ERROR COMPRESS!!!!!!!\n";
  //       }
  //       // if(((*nb)!=(it->first.first)) || ((*pd)!=(it->first.second))) std::cout << "ERROR COMPRESS!!!!!!!\n";
  //     }
  //     // std::cout<< "===============================================\n";
  //   }
  // }

  // for(int i=0; i<new_size; ++i){
  //   for(auto it: cnb[i]){
  //     cnb2[origianl_ref[i]].push_back(origianl_ref[it]);
  //   }
  // }
  
  // printf("|E| after compressing is: %u\n", cedges);
  // printf("|E| after compressing is (indeg check): %u\n", cedges2);
}

void Index::outputGraph(){
  std::ofstream outfile("./tolc/twitter_final_c.in");
  // std::ofstream outfile("./compressed graph/citeseerx_tol_original_3.in");

  if (!outfile.is_open()) {
      std::cerr << "Failed to open file for writing." << std::endl;
  }


  //! shuffle G^c
  // outfile << cnb.size()<< " " << "1111";
  // outfile << cnb.size();
  // outfile << std::endl;
  // for (size_t i = 0; i < cnb.size(); ++i) {
  //     outfile << i << ": ";//no need for PLL
  //     // if(i>=cnb.size()){
  //     //   outfile << "-1";
  //     //   outfile << std::endl;
  //     //   continue;
  //     // }
  //     for (size_t j = 0; j < cnb[i].size(); ++j) {
  //         outfile << cnb[i][j];
  //         if (j < cnb[i].size() - 1) {
  //             outfile << " ";
  //         }
  //     }

  //     if(cnb.size()!=0) outfile << " ";
  //     outfile << "-1";

  //     outfile << std::endl;
  // }

  //! G delete style
  // outfile << n_<< " " << "1";
  // outfile << std::endl;
  // for (size_t i = 0; i < n_; ++i) {
  //     outfile << i << " ";
  //     if(redundant[i]) {
  //       outfile << std::endl;
  //       continue;
  //     }
  //     else{
  //       std::vector<unsigned> * nb = g->get_neighbors(i);
  //       for (size_t j = 0; j < nb->size(); ++j) {
  //           if(redundant[(*nb)[j]]==1) continue;
  //           outfile << (*nb)[j];
  //           if (j < nb->size() - 1) {
  //               outfile << " ";
  //           }
  //       }
  //       outfile << std::endl;
  //     }
  // }
  outfile.close();
}

void Index::outputQuery(std::vector<std::pair<unsigned, unsigned> > * queries){
  std::ofstream outfile("./query/query.txt");
  std::ofstream outfilec("./query/query_c.txt");
  if (!outfilec.is_open()) {
      std::cerr << "Failed to open file for writing." << std::endl;
  }

  for (size_t i = 0; i < queries->size(); ++i) {
      // outfile << (*queries)[i].first << " " << (*queries)[i].second << std::endl;
      // outfilec << reduction_ref[(*queries)[i].first] << " " << reduction_ref[(*queries)[i].second] << std::endl;
      outfile << levelNodes[1][(*queries)[i].first%levelNodes[1].size()] << " " << levelNodes[2][(*queries)[i].second%levelNodes[2].size()] << std::endl;
      outfilec << reduction_ref[levelNodes[1][(*queries)[i].first%levelNodes[1].size()]] << " " << reduction_ref[levelNodes[2][(*queries)[i].second%levelNodes[2].size()]] << std::endl;
  }
  outfile.close();
  outfilec.close();
}
//!indexing on the compressed graph

void Index::calculate_order_c(){

  // std::vector<unsigned> * indegree = &cin_degree;
  // std::vector<unsigned> * degree = &cout_degree;
  int outdeg_ex, indeg_ex;

  for (unsigned i=0; i<cin_degree.size(); ++i){
    outdeg_ex=0;
    indeg_ex=0;
    for(auto & it: cnb[i]) outdeg_ex += ref_clusters[it].size();
    for(auto & it: cpd[i]) indeg_ex += ref_clusters[it].size();

    double temp_score = (double)(outdeg_ex) * (double)(indeg_ex) / ((double)(outdeg_ex) + (double)(indeg_ex) + 1);//new
    if (temp_score < 0) std::cout<< "OVERFLOW" << std::endl;
    
    queue_unit temp(i, temp_score, outdeg_ex+indeg_ex);
    pq_node_score.push(temp);
  }

}

void Index::_2hop_c(){

  std::vector<unsigned> * indegree = &cin_degree;
  std::vector<unsigned> * degree = &cout_degree;

  std::vector<unsigned> q;

  visited = std::vector<int>(CIN.size(),0);

  int i = 0;
  unsigned currentNode;

  while (!pq_node_score.empty()){

    // Timer t; t.start();
    while(1){
      int top_node = pq_node_score.top().node;
      double top_score = pq_node_score.top().score;
      // pq_node_score.pop();
      
      double temp_score = (double)((*indegree)[top_node]) * (double)((*degree)[top_node]) / (double)((*indegree)[top_node] + (*degree)[top_node] + 1) ;//new

      if (temp_score < 0) std::cout<< "OVERFLOW" << std::endl;

      // std::cout<< top_node << "; " << temp_score << "; " << top_score << std::endl;
      if (temp_score == top_score){
        currentNode = top_node;
        pq_node_score.pop();
        break;
      }
      pq_node_score.pop();
      queue_unit temp_struct(top_node, temp_score, (*indegree)[top_node] + (*degree)[top_node]);
      pq_node_score.push(temp_struct);
    }
    // debug_time += t.stop();
    // std::cout<< "choosen node: " << currentNode << "; indeg: " << (*indegree)[currentNode] << "; outdeg: " << (*degree)[currentNode] << std::endl;
    // std::cout<< "level: " << level[currentNode] << std::endl;//error information sincce currentnode changes.
    if ((*indegree)[currentNode] + (*degree)[currentNode] == 0) continue;
    
    ++i;
  
    // Timer t; t.start();
    q.clear();
    q.push_back(currentNode);

    for (unsigned z=0; z<q.size(); ++z){
      unsigned temp = q[z];

      if (visited[temp]==i || visited[temp]==-1) continue;//visited before

      visited[temp]=i;

      if (temp==currentNode) {
        CIN[temp].push_back(i);
        // labelsize += 1;

        std::vector<unsigned> *nb = &cnb[temp];

        for (unsigned k=0; k<nb->size(); ++k){
          q.push_back((*nb)[k]);
        }

      }

      else if (!mergeCheck_c(&currentNode, &temp)) {
        CIN[temp].push_back(i);
        labelsize_c += 1;

        std::vector<unsigned> *nb = &cnb[temp];

        for (unsigned k=0; k<nb->size(); ++k){
          q.push_back((*nb)[k]);
        }

      }
    }
    // debug_time += t.stop();
    // Timer t; t.start();
    visited[currentNode]=0;
    q.clear();
    q.push_back(currentNode);

    for (unsigned z=0; z<q.size(); ++z){
      unsigned temp = q[z];

      if (visited[temp]==i || visited[temp]==-1) continue;//visited before

      visited[temp]=i;

      if (temp == currentNode) {
        COUT[temp].push_back(i);
        // labelsize += 1;
        std::vector<unsigned> *pd = &cpd[temp];

        for (unsigned k=0; k<pd->size(); ++k){
          q.push_back((*pd)[k]);
        }

      }

      else if (!mergeCheck_c(&temp, &currentNode)) {
        COUT[temp].push_back(i);
        labelsize_c +=1;
        
        std::vector<unsigned> *pd = &cpd[temp];

        for (unsigned k=0; k<pd->size(); ++k){
          q.push_back((*pd)[k]);
        }

      }
    }
    // debug_time += t.stop();

    visited[currentNode] = -1;//set to 2 to permanently mask.

    //4. Delete the information of current node. manipulate on the original graph
    const std::vector<unsigned> * nb = &cnb[currentNode];
    for (unsigned h=0; h<nb->size(); ++h){
      if ((*indegree)[(*nb)[h]] > 0) (*indegree)[(*nb)[h]] -= 1;
    }  
    

    const std::vector<unsigned> * pd = &cpd[currentNode];
    for (unsigned h=0; h<pd->size(); ++h){
      (*degree)[(*pd)[h]] -= 1;
    }
    
  }
  reachable=0;
}

void Index::masknodes(int masksize){
  // std::vector<unsigned> * nb;
  // std::vector<unsigned> * pd;
  visited = std::vector<int>(n_,0);
  int random_integer;
  int lowest=0, highest=n_-1;
  int range=(highest-lowest)+1;
  std::vector<unsigned> choosen(n_,0);
  for(int i=0; i<masksize; ++i){
    while(1){
      random_integer = lowest + rand() % range;
      if(choosen[random_integer]==0) break;
    }
    
    delete_order.push_back(random_integer);
    choosen[random_integer]=1;
  }

  // output updates in dagger format;
  // std::ofstream outfile("./updates/test_origin_dagger.up");

  // if (!outfile.is_open()) {
  //     std::cerr << "Failed to open file for writing." << std::endl;
  // }
  // std::unordered_set<unsigned> temp;
  // for (size_t i = 0; i < delete_order.size(); ++i) {
  //     outfile << "r " << delete_order[i] << std::endl;
  //     temp.insert(delete_order[i]);
  // }
  // for (int i = delete_order.size()-1; i>=0; --i) {
  //     nb = g->get_neighbors(delete_order[i]);
  //     pd = g->get_predecessors(delete_order[i]);
  //     temp.erase(delete_order[i]);
  //     int insize = pd->size();
  //     for(auto it: *pd){
  //       if(temp.count(it)) insize--;
  //     }
  //     outfile << "a " << delete_order[i];
  //     outfile << " " << insize << " ";
  //     for(auto it: *pd){
  //       if(!temp.count(it)) outfile << it << " ";
  //     }
  //     int outsize = nb->size();
  //     for(auto it: *nb){
  //       if(temp.count(it)) outsize--;
  //     }
  //     outfile << outsize << " ";
  //     for(auto it: *nb){
  //       if(!temp.count(it)) outfile << it << " ";
  //     }
  //     outfile << std::endl;
  // }
}

void Index::vector_delete(std::vector<unsigned> * vec, unsigned target){
  std::vector<unsigned>::iterator position = std::lower_bound(vec->begin(), vec->end(), target);
  if(*position!= target) std::cout << *position << ";" << target << "Binary Search Error...\n";
  else vec->erase(position);
}

void Index::vector_insert(std::vector<unsigned> * vec, unsigned target){
  std::vector<unsigned>::iterator position = std::lower_bound(vec->begin(), vec->end(), target);
  if(*position!= target) std::cout << *position << ";" << target << "Binary Search Error...\n";
  else vec->erase(position);
}

void Index::maintain(){
  //lazy update in the end: recompress and record the first ref.
  int save = 0;
  // int mergetime = 0;
  // std::vector<unsigned> nodes_to_delete;
  std::vector<unsigned> * nb;
  std::vector<unsigned> * pd;
  
  std::vector<unsigned> * nb_update;
  std::vector<unsigned> * pd_update;
  // int i=0;
  std::unordered_set<unsigned> delete_nodes;
  for(auto node: delete_order){
    // std::cout << i++ << std::endl;
    delete_nodes.insert(node);
    nb = g->get_neighbors(node);
    pd = g->get_predecessors(node);
    // vectors = std::make_pair(*nb, *pd);
    if(ref_clusters[reduction_ref[node]].size() > 1){
      save++;
      ref_clusters[reduction_ref[node]].erase(node);
    }
    else{//==,the only one, delete corresp node.
      ref_clusters[reduction_ref[node]].erase(node);
      // outfile << reduction_ref[node] << std::endl;
    }
    //now update neighnours for merging
 
    for(auto pdnode: *pd){
      visited[pdnode]=1;
      nb_update = g->get_neighbors(pdnode);
      vector_delete(nb_update, node);
    }
    for(auto nbnode: *nb){
      visited[nbnode]=1;
      pd_update = g->get_predecessors(nbnode);
      vector_delete(pd_update, node);
    }
  }
  //lazy update
  for(unsigned i=0; i<n_; i++){
    if(!delete_nodes.count(i) && visited[i]==1){
      // std::cout<< c++ << std::endl;
      std::vector<unsigned> * nb = g->get_neighbors(i);
      std::vector<unsigned> * pd = g->get_predecessors(i);
      vectors = std::make_pair(*nb, *pd);
      std::unordered_set<unsigned> * temp = &(vectorGroups[vectors]);
      if(temp->size()>0){
        reduction_ref[i] = reduction_ref[*(temp->begin())];
        temp->insert(i);
      } 
      else{
        temp->insert(i);
      }
    }
  }
  std::cout<< "Save " << save << " of total " << delete_order.size() << " nodes" << std::endl;
}