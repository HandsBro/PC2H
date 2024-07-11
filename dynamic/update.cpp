#include <cstdio>
#include "utility1.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>
using namespace std;


const int VectorDefaultSize = 20;
const int TOPNUM = 1;

iVector<int> query_ans;

struct Butterfly
{
	int N, lcnt, dir, OPR;

	iVector<int> q,level,list[2],l2n,n2x[2],rnodes,order[2];

	iVector<int> nodes;

	iVector<double> cost[2];

	iVector<Triple<int> > candi;

	iMap<int> mark,cleanmark,random_nodes;

	iHeap<double> pq;

	iHeap<int> opq;

	clock_t tstart,tstop;

	int* sum;

	iVector<iVector<int> > labels[2], backlabels[2], links[2];

	iVector<bool> modified;

	LevelHash levelhash;

        double empty, nonempty, bnonempty, bempty;

	//iMap<int> father;

	inline int AddLabel( int s, int t, int d, int index )
	{
//                if ( labels[d][s].m_num > 100 && labels[1-d][t].m_num > 100 ) 
//                   printf("%d(%d) %d(%d)\n", s, labels[d][s].m_num, t, labels[1-d][t].m_num);
//
//              

/*                bool flag = true;
                for ( int i = 0 ; i < 5 ; ++i )
                {                
                    if ( (labels[d][s].mask[i] & labels[1-d][t].mask[i]) != 0 )
                    {
                            flag = false;
                            break;
                   }
                }

                if ( flag )
                {
                            labels[d][s].push_back_with_mask( index );
                            return 1;                
                }
           
                int ns = labels[d][s].m_num, nt = labels[1-d][t].m_num;

                if ( ns > 0 && nt > 0 && ns + nt > 1000 )
                {
                        float c1 = log((float)ns)*10.0*nt, c2 = log((float)nt)*10.0*ns, c3 = ns+nt;

                        if ( c1 <= c2 && c1 <= c3 )
                        {
                            for ( int j = 0 ; j < labels[1-d][t].m_num ; ++j )
                            {
                                if ( labels[d][s].BinarySearch(labels[1-d][t][j]) >= 0 )
                                {
                                        bnonempty += log((float)ns)*10.0*(j+1);
                                        return 0;                            
                                }
                            }

                            labels[d][s].push_back_with_mask(index);
                            bempty += c1;
                            return 1;
                        }
                        

                        if ( c2 <= c1 && c2 <= c3 )
                        {
                            for ( int i = 0 ; i < labels[d][s].m_num ; ++i )
                            {
                                if ( labels[1-d][t].BinarySearch(labels[d][s][i]) >= 0 )
                                {
                                    bnonempty += log((float)nt)*10.0*(i+1);
                                    return 0;                            \
                                }
                            }

                            labels[d][s].push_back_with_mask(index);
                            bempty += c2;
                            return 1;                        
                        }
                        
                }
            
		int i , j , x , y;
		for ( i = 0 , j = 0 ; i < labels[d][s].m_num && j < labels[1-d][t].m_num ; )
		{
			x = labels[d][s][i]; y = labels[1-d][t][j];
			if ( x == y )
                        {
                                nonempty+=(i+j);
                                return 0;
                        }
			else if ( x < y ) i++;
			else j++;
		}
*/
		labels[d][s].push_back( index );
//                empty+=(ns+nt);
		return 1;
	}

	//! not used concrete nodes as lables. Use index instead and it is natually orderred.
	inline int AddLabelX( int s, int t, int d, int index )
	{
		int l = 0 , r = labels[d][s].m_num, o = cleanmark.occur.m_num, p, m, x;

		if ( o > 0 && r > 0 )
		{
			x = r/(o*8);
			if ( x > 30 || r < (1<<(x)) )
			{
				for ( int j = 0 ; j < o ; r = labels[d][s].m_num, j++ )
				{
					p = cleanmark.occur[j];
					for ( ; l < r ; )
					{
						m = (l+r)/2;
						x = labels[d][s][m];

						if ( x == p ) return 0;
						else if ( x < p ) l = m+1;
						else r = m;
					}
				}
			}
			else
			{
				for ( int i = 0 ; i < r ; ++i )
					if ( cleanmark.exist(labels[d][s][i]) ) return 0;
			}
		}

		labels[d][s].push_back( index );
		return 1;
	}


	inline int AddLabel1( int s, int t, int d )
	{
		int i , j , x , y;
		for ( i = 0 , j = 0 ; i < labels[d][s].m_num && j < labels[1-d][t].m_num ; )
		{
			x = labels[d][s][i]; y = labels[1-d][t][j];
			if ( x == y ) return 0;
			else if ( x < y ) i++;
			else j++;
		}

		labels[d][s].sorted_insert(t);

		backlabels[d][t].sorted_insert(s);
		
		////insert label from s to t
		//labels[d][s].push_back( t );
		//for ( int i = labels[d][s].m_num-1 ; i > 0 ; --i )
		//{
		//	if ( labels[d][s][i] < labels[d][s][i-1] )
		//	{
		//		int tmp = labels[d][s][i];
		//		labels[d][s][i] = labels[d][s][i-1];
		//		labels[d][s][i-1] = tmp;
		//	}
		//	else break;
		//}
		////insert backlabel from t to s
		//backlabels[d][t].push_back(s);
		//for ( int i = backlabels[d][t].m_num-1 ; i > 0 ; --i )
		//{
		//	if ( backlabels[d][t][i] < backlabels[d][t][i-1] )
		//	{
		//		int tmp = backlabels[d][t][i];
		//		backlabels[d][t][i] = backlabels[d][t][i-1];
		//		backlabels[d][t][i-1] = tmp;
		//	}
		//	else break;
		//}		
		return 1;
	}

	void query( string query_file, int code )
	{
		const int testcases = 1000000;

		int *queries = new int[testcases*2];		

		FILE *file = fopen( query_file.c_str(), "r" );

		for ( int i = 0 ; i < testcases ; ++i )
		{
			fscanf( file, "%d %d", &queries[i*2], &queries[i*2+1] );
		}

		fclose(file);

		int sum = 0, s, t, x, y, i, j;

	//	Runtimecounter rt;

	//	rt.start();
		
		for ( int p = 0 ; p < testcases*2 ; ++p )
		{
			s = queries[p++]; t = queries[p];

			//if ( p/2 == 605337 )
			//{
			//	printf("shit %d %d\n",s,t);

			//	if ( random_nodes.occur.m_num > 0 )
			//	{
			//		for ( int i = 0 ; i < labels[0][s].m_num ; ++i )
			//		{
			//			printf("ooo %d(%d) %d(%d)\n",s,random_nodes.get(s),labels[0][s][i],random_nodes.get(labels[0][s][i]));
			//		}
			//		for ( int i = 0 ; i < labels[1][t].m_num ; ++i )
			//		{
			//			printf("ooo %d(%d) %d(%d)\n",t,random_nodes.get(t),labels[1][t][i],random_nodes.get(labels[1][t][i]));
			//		}
			//	}
			//	else
			//	{
			//		for ( int i = 0 ; i < labels[0][s].m_num ; ++i )
			//		{
			//			printf("nnn %d %d\n",s,l2n[labels[0][s][i]]);
			//		}
			//		for ( int i = 0 ; i < labels[1][t].m_num ; ++i )
			//		{
			//			printf("nnn %d %d\n",t,l2n[labels[1][t][i]]);
			//		}
			//	}
			//}

			//if ( p/2 == 546808 )
			//{
			//	printf("shit %d %d\n",s,t);

			//	if ( random_nodes.occur.m_num > 0 )
			//	{
			//		for ( int i = 0 ; i < labels[0][s].m_num ; ++i )
			//		{
			//			printf("new %d %d\n",s,labels[0][s][i]);
			//		}
			//		for ( int i = 0 ; i < labels[1][t].m_num ; ++i )
			//		{
			//			printf("new %d %d\n",t,labels[1][t][i]);
			//		}
			//	}
			//	else
			//	{
			//		for ( int i = 0 ; i < labels[0][s].m_num ; ++i )
			//		{
			//			printf("old %d %d\n",s,l2n[labels[0][s][i]]);
			//		}
			//		for ( int i = 0 ; i < labels[1][t].m_num ; ++i )
			//		{
			//			printf("old %d %d\n",t,l2n[labels[1][t][i]]);
			//		}
			//	}
			//}

			for ( i = 0 , j = 0 ; i < labels[0][s].m_num && j < labels[1][t].m_num ; )
			{
				x = labels[0][s][i]; y = labels[1][t][j];
				if ( x == y ) 
				{
					sum++;
/*					if ( code == 1 ) query_ans.push_back( 1 );
					else if ( code == 2 && query_ans[p/2] != 1 )
					{
						if ( labels[0][s].m_num + labels[1][t].m_num > 0 ) 
						{
							printf("Error1: %d %d %d %d %d\n",p/2,s,t,labels[0][s].m_num,labels[1][t].m_num);
							return;
						}
					}*/
					break;
				}
				else if ( x < y ) i++;
				else j++;
			}

/*			if ( code == 1 && query_ans.m_num <= p/2 ) query_ans.push_back(0);
			else if ( code == 2 &&  ( i == labels[0][s].m_num || j == labels[1][t].m_num ) && query_ans[p/2] == 1 )
			{
						if ( labels[0][s].m_num + labels[1][t].m_num > 0 ) 
						{
							printf("Error2: %d %d %d %d %d\n",p/2,s,t,labels[0][s].m_num,labels[1][t].m_num);
							//return;
						}
			}*/
		}

//		rt.stop();

		printf("%d\n",sum);

//		printf("Time: %0.8lf\n",rt.GetRuntime());
	}

	bool circle_check( int s, int t, int d )
	{
		for ( int i = 0 ; i < labels[1-d][t].m_num ; ++i )
		{
			if ( labels[1-d][t][i] == s )
			{
				int x , y, cnt = 0;
				for ( int j = 0 , k = 0 ; j < backlabels[1-d][s].m_num && k < backlabels[1-d][t].m_num ; )
				{
					x = backlabels[1-d][s][j]; y = backlabels[1-d][t][k];
					if ( x == y )
					{
						cnt++;
						++k;
						++j;
					}
					else if ( x < y ) j++;
					else k++;
				}
				if ( cnt == backlabels[1-d][t].m_num ) 
				{
					return true;
				}
				else return false;
			}
		}
		return false;
	}

	bool bicheck( int p, int d, int t )
	{
		int m, l = 0, r = labels[d][t].m_num;
		for ( ; l < r ; )
		{
			m = (l+r)/2;
			if ( labels[d][t][m] == p ) return true;
			if ( labels[d][t][m] < p ) l = m+1;
			else r = m;
		}
		return false;
	}

	void load_tflabel( string st )
	{
		int x, y = 0;

		//load level
		string levelfn = "index\\"+st+"_level";
		FILE *Lfile = fopen(levelfn.c_str(),"r");

		mark.clean();

		fscanf(Lfile,"%d",&x);

		for ( int i = 0 ; i < x ; ++i )
		{

			fscanf(Lfile,"%d",&y);
			level[y] = i;
			l2n[i] = y;
			mark.insert(y,-2);
		}

		for ( int i = 0 ; i < N ; ++i )
		{
			if ( fscanf(Lfile,"%d",&y) != 1 ) printf("Error!!!\n");
			if ( mark.notexist(i) ) mark.insert(i,y);
		}

		fclose(Lfile);

		y = N;
		
		for ( int r = 0 ; r < N ; ++r )
		{
			bool flag = true;

			for ( int i = N-1 ; i >= 0 ; --i )
			{
				if ( mark.get(i) == r )
				{
					--y;
					flag = false;
					level[i] = y;
					l2n[y] = i;
				}
			}

			if ( flag )
			{
				for ( int i = N-1 ; i >= 0 ; --i )
				{
					if ( mark.get(i) == -1 )
					{
						--y;
						level[i] = y;
						l2n[y] = i;
					}
				}
				break;
			}
		}

		if ( x != y ) printf("Error!!!\n");


		//load labels
		iVector<int> start;
		string startfn = "index\\"+st+"_tlstart";
		FILE* file = fopen(startfn.c_str(),"rb");
		for ( y = 0 ; fread(&x,sizeof(int),1,file) == 1 ; )
		{
			y += x;
			start.push_back(y);
		}
		fclose(file);

		if ( N != start.m_num/2 ) printf("Something Weird!!!\n");

		string labelfn = "index\\"+st+"_TL";
		FILE* lfile = fopen( labelfn.c_str(), "rb" );

		lcnt = 0;

		backlabels[0].re_allocate(N); backlabels[0].m_num = N;
		backlabels[1].re_allocate(N); backlabels[1].m_num = N;

		int errorcnt = 0;

		for ( int i = 0 , p = 0 ;  fread(&x,sizeof(int),1,lfile) == 1 ; ++i )
		{
			if ( i == start[p] ) p++;			
			if ( level[x] <= level[p%N] )
			{
				labels[p/N][p%N].push_back(x);
				backlabels[p/N][x].push_back(p%N);
				if ( x != p%N ) 
				{
					lcnt++;
				}
			}
			else
			{
				errorcnt++;
				if ( p%N == 333945 && x == 18186 )
				{
					printf("fuck\n");
				}
			}
		}

		printf("%d correct TF labels loaded, %d wrong labels, start computing levels...\n",lcnt,errorcnt);


		//Triangle Reduce
		mark.clean();
		
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			n2x[dir].re_allocate(N); n2x[dir].m_num = N;
			for ( int i = 0 ; i < N ; ++i )
			{
				n2x[dir][i] = (int)(log((double)(labels[dir][i].m_num))/log(2.0)+1);
			}
		}

		int coef = 84;
		iVector<int> tmp;

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < N ; ++i )
			{
				if ( labels[dir][i].m_num > 1 )
				{
					//if ( start[dir][i+1]-start[dir][i] > 10000 ) printf("%d %d %d\n",dir,i,start[dir][i+1]-start[dir][i]);
					for ( int j = 0 ; j < labels[dir][i].m_num ; ++j )
					{
						if ( labels[dir][i][j] != i ) 
						{
							mark.insert(labels[dir][i][j],1);
						}
					}
					for ( int j = 0 ; j < labels[dir][i].m_num ; ++j )
					{
						int t = labels[dir][i][j];
						bool flag = true;
						if ( t != i )
						{
							if ( labels[1-dir][t].m_num > labels[dir][i].m_num*coef*n2x[1-dir][t] )
							{
								for ( int h = 0 ; h < labels[dir][i].m_num ; ++h )
								{
									if ( labels[dir][i][h] != i && labels[dir][i][h] != t && bicheck( labels[dir][i][h], 1-dir, t ) )
									{
										flag = false;
										break;
									}
								}
							}
							else
							{
								for ( int k = 0 ; k < labels[1-dir][t].m_num ; ++k )
								{
									if ( labels[1-dir][t][k] != t && mark.exist(labels[1-dir][t][k]) )
									{
										flag = false;
										break;
									}
								}
							}
						}
						if ( flag ) tmp.push_back(t);
						else
						{
							if ( !backlabels[dir][t].remove( i ) ) 
							{
								printf("Backlabel deleting Error!!!\n");
							}
							lcnt--;
						}
					}
					mark.clean();
					labels[dir][i].clean();
					labels[dir][i].push_back( tmp.m_data, tmp.m_num );
					tmp.clean();
				}
			}
		}

		printf("%d labels after Triangle Reduce\n",lcnt);
	}

	void optimize()
	{
		levelhash.initialize( level.m_data, N, 30 );

		tstart = clock();

		printf("\n");

		int sum = 0;
		OPR = 0;
		for ( int i = 0 ; i < N-TOPNUM ; ++i )
		{
			if ( i % 1000 == 0 )
			{
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d: %0.0lf + %0.0lf = %0.0lf",i,4.0*lcnt,4.0*sum,4.0*(lcnt+sum));
			}
			sum += UpgradeNode( l2n[N-1-i] );
		}		

		tstop = clock();

		printf( "\nOptimizing operations: %d with Time: %0.2lf\n",OPR,(double)(tstop-tstart)/CLOCKS_PER_SEC );

		printf("Reconstruct Times: %d\n", levelhash.R);

	}

	void load_graph( string graph_name, int extra )
	{
		FILE *file = fopen(graph_name.c_str(),"r");

		int useless = fscanf(file,"%d",&N);

		printf("%d nodes\n",N);

		mark.initialize(N+extra);
		cleanmark.initialize(N+extra);
		pq.initialize(N);

		for ( dir = 0 ; dir < 2 ; ++dir )//0:nb; 1:pd.
		{
			labels[dir].re_allocate(N+extra);
			labels[dir].m_num = N;//why set to the last with empty bucket before?: here refers to the outter container. 
			cost[dir].re_allocate(N);	
			cost[dir].m_num = N;
			links[dir].re_allocate(N+extra);
			links[dir].m_num = N;
		}
		q.re_allocate(N+extra);
		level.re_allocate(N+extra); level.m_num = N+extra;
		l2n.re_allocate(N+extra); l2n.m_num = N+extra;

		sum = new int[N+1];
		memset(sum,0,(N+1)*sizeof(int));

		int x, y, lines = 0;
		for ( ; fscanf(file,"%d:",&x) == 1 ; ++lines)
		{		
			for ( ; fscanf(file,"%d",&y) == 1 ; )
			{
				// cout<< x << "; " << y << endl;
				if ( y == -1 ) break;
				int t = links[1][y].m_num, c = links[0][x].m_num;
				// cout << t << "; " << c << endl;
				links[0][x].push_back(y);//nb
				links[1][y].push_back(x);//pd
			}
		}

		printf("%d lines scaned\n",lines);

		fclose(file);
	}

	void compute_index_r1()
	{
		// compute_order();
		cleanmark.initialize(N);
		double costlimit = 1e20;//double discard the overflow problem.
        //empty = nonempty = bempty = bnonempty = 0;
               
		tstart = clock();

//		father.initialize(N);

		//compute topological order and cost
		//cost is the approximate in/out set size of each node.

		for ( dir = 0 ; dir < 2 ; ++dir )//0: nb. 1: pd. 
		{
			for ( int i = 0 ; i < N ; ++i )
			{
				sum[i] = links[1-dir][i].m_num;//
				if ( sum[i] == 0 ) q.push_back(i);
				cost[dir][i] = 1;//set all to 1 as init for leaves and roots. 
			}		

			for ( int i = 0 ; i < q.m_num ; ++i )//
			{
				int p = q[i];

				for ( int j = 0 ; j < links[dir][p].m_num ; ++j )
				{
					int t = links[dir][p][j];//
					cost[dir][t] += cost[dir][p];//+1
					if ( cost[dir][t] > costlimit )
                                        {
                                                cost[dir][t] = costlimit*costlimit;
												// cout<< cost[dir][t] << "@@@@@@@@@@@@@\n";
                                        }
					sum[t]--;
					if ( sum[t] == 0 ) q.push_back(t);//BFS models.
				}
			}

			q.clean();
		}

		//compute hierarchies
		lcnt = 0;//label counter
		int ins=0;
		int outs=0;

		double *tcost = new double[N];

		for ( int i = 0 ; i < N ; ++i )
		{
			// cout<< i << ": " << cost[0][i] << "; " << cost[1][i] << ". costlimit: " << costlimit <<endl;
                        if ( cost[0][i] > costlimit || cost[1][i] > costlimit )//if too large.
                        {
                            pq.insert(i,-costlimit*costlimit*(links[0][i].m_num+1)*(links[1][i].m_num+1));//negative is used for pq.
                        }
                        else
                        {            
			    			pq.insert(i,-(cost[0][i]-1)*(cost[1][i]-1)/(cost[0][i]+cost[1][i]));
							// cout<< i << "; " << -(cost[0][i]-1)*(cost[1][i]-1)/(cost[0][i]+cost[1][i]) << endl;
                        }
			level[i] = N;
			// cout<< level[i] << endl;
		}

		for ( int id = 0 ; id < N ; ++id )
		{

			// if ( id % 100 == 0 )//== 19195381 ) 
			// if (id==10000)
			// {
			// 	printf("%d\n",id);
			// 	break;
			// }
			//get top node with lazy update
			int p ;

			for ( ; true ; )
			{
				p = pq.head();
                                if ( cost[0][p] > costlimit || cost[1][p] > costlimit )
                                {
                                    pq.insert(p,-costlimit*costlimit*(links[0][p].m_num+1)*(links[1][p].m_num+1));
                                }
                                else pq.insert( p, -(cost[0][p]-1)*(cost[1][p]-1)/(cost[0][p]+cost[1][p]) );
				if ( pq.head() == p )
				{
					p = pq.pop();
					// cout<< "chosen node: " << p << "; indeg: " << links[1][p].m_num  << "; deg: " << links[0][p].m_num << "; orders: " << order[1][p]  << endl;
					break;
				}
			}

		//	father.clean();


			//find children of p, mark them, and create back labels.
			//Back label refers to IN labels?
			// cout<< p << "; costs: " << cost[0][p] << "; " << cost[1][p] << endl;
			q.push_back(p);
			mark.insert(p,2);//what does this mean?: means p is accessed downward.
			tcost[p] = 0;

//                        if ( id % 10000 == 0 ) 
//                            printf("%d %0.2lf %0.2lf %0.2lf %0.2lf\n",id,empty/id,nonempty/id,bempty/id,bnonempty/id);

			cleanmark.clean();//mark the out labels of chosen node p.

			for ( int i = 0 ; i < labels[0][p].m_num ; ++i )
			{
				cleanmark.insert( labels[0][p][i], 1 );
			}

			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];

//                                if ( AddLabel( c, p, 1, id ) && c != p ) lcnt++;
 
				if ( AddLabelX( c, p, 1, id ) )//what does this mean?:whether add into in labels. 1 refers to in labels. 
				{
					// lcnt++;
					if ( c != p ) {
						lcnt++;
						ins++;
					}
					for ( int j = 0 ; j < links[0][c].m_num ; ++j )
					{
						int t = links[0][c][j];
						if ( level[t] == N )
						{
							if ( mark.notexist(t) )//if not visited before.
							{
								mark.insert(t,1);//downward accessed nodes are marked as 1. 
								q.push_back(t);		
								tcost[t] = 0;
							}
						}
					}
				}
			}
			q.clean();



			//find ancestors of p
			cleanmark.clean();

			for ( int i = 0 ; i < labels[1][p].m_num ; ++i )
			{
				cleanmark.insert( labels[1][p][i], 1 );
			}

			q.push_back(p);
			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];

				if ( AddLabelX( c, p, 0, id ) )
				{
					// lcnt++;
					if ( c != p ) {
						lcnt++;//self do not count.
						outs++;
					}
					for ( int j = 0 ; j < links[1][c].m_num ; ++j )
					{
						int t = links[1][c][j];
						if ( level[t] == N )
						{
							if ( mark.notexist(t) )
							{
								mark.insert(t,2);//upward accessed nodes marks as 2: ancestors
								q.push_back(t);
								tcost[t] = 0;
							}
							else
							{
								if ( mark.get(t) == 1 )//implies cycle.
								{
									printf("DAG Property Error!!!\n");
								}
							}
						}
					}
				}
			}
			q.clean();

		
			//if ( id == 0 )
			//{
			//	printf("ffff\n");
			//}

			for ( int i = 0 ; i < mark.occur.m_num ; ++i )
			{
				int c = mark.occur[i];
				if ( mark.get(c) == 2 )//ancestors
				{
					for ( int j = 0 ; j < links[0][c].m_num ; ++j )//nb
					{
						int t = links[0][c][j];
						if ( level[t] == N && ( t == p || mark.get(t)==1 ) )
						{
							tcost[c] += cost[1][t];//upward cost?
							if ( tcost[c] > costlimit ) tcost[c] = costlimit;
							tcost[t] += cost[0][c];//down cost?
							if ( tcost[t] > costlimit ) tcost[t] = costlimit;
						}
					}
				}
			}

			for ( int i = 1 ; i < mark.occur.m_num ; ++i )// since the first one is p itself.
			{
				int c = mark.occur[i];
				int d = mark.get(c)-1;
				if ( cost[d][c] < costlimit ) 
                                    cost[d][c] -= tcost[c];

				for ( int j = 0 ; j < links[d][c].m_num ; ++j )
				{
					int t = links[d][c][j];

					if ( level[t] == N )
					{
						tcost[t] += tcost[c];
						if ( tcost[t] > costlimit ) tcost[t] = costlimit;
					}
				}
			}

			level[p] = id;			

			l2n[id] = p;

			mark.clean();
		}

	//	reduce();

		tstop = clock();
		// printf( "Time: %0.5lf ms\n",(double)(tstop-tstart) * 1000 /CLOCKS_PER_SEC );
		// printf( "Label Size: %0.0lf\n",lcnt*4.0);
		// printf( "IN Size: %d\n",ins);
		// printf( "OUT Size: %d\n",outs);
		delete[] tcost;
	}
       

	void compute_index_p1()
	{
		cleanmark.initialize(N);
		double costlimit = 1e20;
        //empty = nonempty = bempty = bnonempty = 0;
               
		tstart = clock();

//		father.initialize(N);

		//compute topological order and cost
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < N ; ++i )
			{
				sum[i] = links[1-dir][i].m_num;
				if ( sum[i] == 0 ) q.push_back(i);
				cost[dir][i] = 1;
			}		

			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int p = q[i];

				for ( int j = 0 ; j < links[dir][p].m_num ; ++j )
				{
					int t = links[dir][p][j];
					cost[dir][t] += cost[dir][p]/links[dir][p].m_num;
					if ( cost[dir][t] > costlimit )
                                        {
                                                cost[dir][t] = costlimit*costlimit;
                                        }
					sum[t]--;
					if ( sum[t] == 0 ) q.push_back(t);
				}
			}

			q.clean();
		}

		//compute hierarchies
		lcnt = 0;

		double *tcost = new double[N];

		for ( int i = 0 ; i < N ; ++i )
		{
                        if ( cost[0][i] > costlimit || cost[1][i] > costlimit )
                        {
                            pq.insert(i,-costlimit*costlimit*(links[0][i].m_num+1)*(links[1][i].m_num+1));
                        }
                        else
                        {            
			    pq.insert(i,-(cost[0][i]-1)*(cost[1][i]-1)/(cost[0][i]+cost[1][i]));
                        }
			level[i] = N;
		}

		for ( int id = 0 ; id < N ; ++id )
		{

			//if ( id % 1000 == 0 )//== 19195381 ) 
			//{
			//	printf("%d\n",id);
			//}
			//get top node with lazy update
			int p ;

			for ( ; true ; )
			{
				p = pq.head();
                                if ( cost[0][p] > costlimit || cost[1][p] > costlimit )
                                {
                                    pq.insert(p,-costlimit*costlimit*(links[0][p].m_num+1)*(links[1][p].m_num+1));
                                }
                                else pq.insert( p, -(cost[0][p]-1)*(cost[1][p]-1)/(cost[0][p]+cost[1][p]) );
				if ( pq.head() == p )
				{
					p = pq.pop();
					break;
				}
			}

		//	father.clean();

			// cout<< "now at node: " << p << endl;
			//find children of p, mark them, and create back labels
			q.push_back(p);
			mark.insert(p,2);
			tcost[p] = 0;

//                        if ( id % 10000 == 0 ) 
//                            printf("%d %0.2lf %0.2lf %0.2lf %0.2lf\n",id,empty/id,nonempty/id,bempty/id,bnonempty/id);

			cleanmark.clean();

			for ( int i = 0 ; i < labels[0][p].m_num ; ++i )
			{
				cleanmark.insert( labels[0][p][i], 1 );
			}

			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];

//                                if ( AddLabel( c, p, 1, id ) && c != p ) lcnt++;
 
				if ( AddLabelX( c, p, 1, id ) )
				{
					if ( c != p ) lcnt++;
					for ( int j = 0 ; j < links[0][c].m_num ; ++j )
					{
						int t = links[0][c][j];
						if ( level[t] == N )
						{
							if ( mark.notexist(t) )
							{
								mark.insert(t,1);
								q.push_back(t);		
								tcost[t] = 0;
							}
						}
					}
				}
			}
			q.clean();



			//find ancestors of p
			cleanmark.clean();

			for ( int i = 0 ; i < labels[1][p].m_num ; ++i )
			{
				cleanmark.insert( labels[1][p][i], 1 );
			}

			q.push_back(p);
			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];

//                                if ( AddLabel( c, p, 0, id ) && c != p ) lcnt++;


				if ( AddLabelX( c, p, 0, id ) )
				{
					if ( c != p ) lcnt++;
					for ( int j = 0 ; j < links[1][c].m_num ; ++j )
					{
						int t = links[1][c][j];
						if ( level[t] == N )
						{
							if ( mark.notexist(t) )
							{
								mark.insert(t,2);
								q.push_back(t);
								tcost[t] = 0;
							}
							else
							{
								if ( mark.get(t) == 1 )
								{
									printf("DAG Property Error!!!\n");
								}
							}
						}
					}
				}
			}
			q.clean();

		
			//if ( id == 0 )
			//{
			//	printf("ffff\n");
			//}

			for ( int i = 0 ; i < mark.occur.m_num ; ++i )
			{
				int c = mark.occur[i];
				if ( mark.get(c) == 2 )
				{
					for ( int j = 0 ; j < links[0][c].m_num ; ++j )
					{
						int t = links[0][c][j];
						if ( level[t] == N && ( t == p || mark.get(t)==1 ) )
						{
							tcost[c] += cost[1][t]/links[1][t].m_num;
							tcost[t] += cost[0][c]/links[0][c].m_num;
							if ( tcost[c] > costlimit ) tcost[c] = costlimit;
							if ( tcost[t] > costlimit ) tcost[t] = costlimit;
						}
					}
				}
			}

			for ( int i = 1 ; i < mark.occur.m_num ; ++i )
			{
				int c = mark.occur[i];
				int d = mark.get(c)-1;
				if ( cost[d][c] < costlimit ) 
                                    cost[d][c] -= tcost[c];

				for ( int j = 0 ; j < links[d][c].m_num ; ++j )
				{
					int t = links[d][c][j];

					if ( level[t] == N )
					{
						tcost[t] += tcost[c]/links[d][c].m_num;
						if ( tcost[t] > costlimit ) tcost[t] = costlimit;
					}
				}
			}

			level[p] = id;			

			l2n[id] = p;

			mark.clean();
		}

	//	reduce();

		tstop = clock();

		printf( "Time: %0.2lf\n",(double)(tstop-tstart)/CLOCKS_PER_SEC );

		printf( "Label Size: %0.0lf\n",lcnt*4.0);

		delete[] tcost;
	}





	void compute_index_r()
	{
		double costlimit = 1e20;
        empty = nonempty = bempty = bnonempty = 0;
               
		tstart = clock();

//		father.initialize(N);

		//compute topological order and cost
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < N ; ++i )
			{
				sum[i] = links[1-dir][i].m_num;
				if ( sum[i] == 0 ) q.push_back(i);
				cost[dir][i] = 1;
			}		

			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int p = q[i];

				for ( int j = 0 ; j < links[dir][p].m_num ; ++j )
				{
					int t = links[dir][p][j];
					cost[dir][t] += cost[dir][p];
					if ( cost[dir][t] > costlimit )
                                        {
                                                cost[dir][t] = costlimit*costlimit;
                                        }
					sum[t]--;
					if ( sum[t] == 0 ) q.push_back(t);
				}
			}

			q.clean();
		}

		//compute hierarchies
		lcnt = 0;

		double *tcost = new double[N];

		for ( int i = 0 ; i < N ; ++i )
		{
                        if ( cost[0][i] > costlimit || cost[1][i] > costlimit )
                        {
                            pq.insert(i,-costlimit*costlimit*(links[0][i].m_num+1)*(links[1][i].m_num+1));
                        }
                        else
                        {            
			    pq.insert(i,-(cost[0][i]-1)*(cost[1][i]-1)/(cost[0][i]+cost[1][i]));
                        }
			level[i] = N;
		}

		for ( int id = 0 ; id < N ; ++id )
		{
			//get top node with lazy update
			int p ;

			for ( ; true ; )
			{
				p = pq.head();
                                if ( cost[0][p] > costlimit || cost[1][p] > costlimit )
                                {
                                    pq.insert(p,-costlimit*costlimit*(links[0][p].m_num+1)*(links[1][p].m_num+1));
                                }
                                else pq.insert( p, -(cost[0][p]-1)*(cost[1][p]-1)/(cost[0][p]+cost[1][p]) );
				if ( pq.head() == p )
				{
					p = pq.pop();
					break;
				}
			}

		//	father.clean();

			//find children of p, mark them, and create back labels
			q.push_back(p);
			mark.insert(p,2);
			tcost[p] = 0;

//                        if ( id % 10000 == 0 ) 
//                            printf("%d %0.2lf %0.2lf %0.2lf %0.2lf\n",id,empty/id,nonempty/id,bempty/id,bnonempty/id);

			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];

//                                if ( AddLabel( c, p, 1, id ) && c != p ) lcnt++;
 

				if ( AddLabel( c, p, 1, id ) )
				{
					if ( c != p ) lcnt++;
					for ( int j = 0 ; j < links[0][c].m_num ; ++j )
					{
						int t = links[0][c][j];
						if ( level[t] == N )
						{
							if ( mark.notexist(t) )
							{
								mark.insert(t,1);
								q.push_back(t);		
								tcost[t] = 0;
							}
						}
					}
				}
			}
			q.clean();

			//find ancestors of p
			q.push_back(p);
			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];

//                                if ( AddLabel( c, p, 0, id ) && c != p ) lcnt++;

				if ( AddLabel( c, p, 0, id ) )
				{
					if ( c != p ) lcnt++;
					for ( int j = 0 ; j < links[1][c].m_num ; ++j )
					{
						int t = links[1][c][j];
						if ( level[t] == N )
						{
							if ( mark.notexist(t) )
							{
								mark.insert(t,2);
								q.push_back(t);
								tcost[t] = 0;
							}
							else
							{
								if ( mark.get(t) == 1 )
								{
									printf("DAG Property Error!!!\n");
								}
							}
						}
					}
				}
			}
			q.clean();

		
			//if ( id == 0 )
			//{
			//	printf("ffff\n");
			//}

			for ( int i = 0 ; i < mark.occur.m_num ; ++i )
			{
				int c = mark.occur[i];
				if ( mark.get(c) == 2 )
				{
					for ( int j = 0 ; j < links[0][c].m_num ; ++j )
					{
						int t = links[0][c][j];
						if ( level[t] == N && ( t == p || mark.get(t)==1 ) )
						{
							tcost[c] += cost[1][t];
							if ( tcost[c] > costlimit ) tcost[c] = costlimit;
							tcost[t] += cost[0][c];
							if ( tcost[t] > costlimit ) tcost[t] = costlimit;
						}
					}
				}
			}

			for ( int i = 1 ; i < mark.occur.m_num ; ++i )
			{
				int c = mark.occur[i];
				int d = mark.get(c)-1;
				if ( cost[d][c] < costlimit ) 
                                    cost[d][c] -= tcost[c];

				for ( int j = 0 ; j < links[d][c].m_num ; ++j )
				{
					int t = links[d][c][j];

					if ( level[t] == N )
					{
						tcost[t] += tcost[c];
						if ( tcost[t] > costlimit ) tcost[t] = costlimit;
					}
				}
			}

			level[p] = id;			

			l2n[id] = p;

			mark.clean();
		}

		reduce();

		tstop = clock();

		printf( "Time: %0.2lf\n",(double)(tstop-tstart)/CLOCKS_PER_SEC );

		printf( "Label Size: %0.0lf\n",lcnt*4.0);

		delete[] tcost;
	}



	void output( string st )
	{
		compute_order();//compute the topological order.

		string s1 = st + ".label";
		string s2 = st + ".index";
		FILE *file = fopen(s1.c_str(),"w");//label
		FILE *ifile = fopen(s2.c_str(),"w");//index
		int maxout=0, maxin=0, outid, inid;


		int x = 0;
		for ( int i = 0 ; i < N ; ++i )
		{
			fprintf(ifile,"%d ",x);

			fprintf(file,"%d ",order[0][i]); 
			x++;
			
			// maxout = max(maxout, (int)labels[0][i].m_num);
			if((int)labels[0][i].m_num>maxout){
				maxout = (int)labels[0][i].m_num;
				outid = i;
			}

			for ( int j = 0 ; j < labels[0][i].m_num ; ++j )
			{
				fprintf(file,"%d ",labels[0][i][j]);
				x++;
			}

			fprintf(file,"%d ",-1);
			x++;

			fprintf(ifile,"%d\n",x);

			fprintf(file,"%d ",order[0][i]); 
			x++;

			// maxin = max(maxin, (int)labels[1][i].m_num);
			if((int)labels[1][i].m_num>maxin){
				maxin = (int)labels[1][i].m_num;
				inid = i;
			}


			for ( int j = 0 ; j < labels[1][i].m_num ; ++j )
			{
				fprintf(file,"%d ",labels[1][i][j]);
				x++;
			}

			fprintf(file,"%d\n",-2);
			x++;
		}

		fclose(file);
		fclose(ifile);
		printf("MAX OUT SIZE: %d; id: %d\n", maxout,outid);
		printf("MAX IN SIZE: %d; id: %d\n", maxin,inid);
	}

	void compute_index_d()
	{
		tstart = clock();

		iVector<Key_Value<double,int> > pn;

		//compute hierarchies
		lcnt = 0;

		for ( int i = 0 ; i < N ; ++i )
		{
			pn.push_back( Key_Value<double,int>( -1.0*(links[0][i].m_num+1)*(links[1][i].m_num+1), i ) );
			level[i] = N;
		}

		pn.Sort();

		for ( int id = 0 ; id < N ; ++id )
		{
			int p = pn[id].value;			

			//find children of p, mark them, and create back labels
			q.push_back(p);
			mark.insert(p,2);

			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];
				if ( AddLabel(c,p,1,id) )
				{
					if ( c!= p ) lcnt++;
					for ( int j = 0 ; j < links[0][c].m_num ; ++j )
					{
						int t = links[0][c][j];
						if ( level[t] == N )
						{
							if ( mark.notexist(t) )
							{
								mark.insert(t,1);
								q.push_back(t);		
							}
						}
					}
				}
			}
			q.clean();

			//find ancestors of p
			q.push_back(p);
			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];
				if ( AddLabel(c,p,0,id) )
				{
					if ( c!= p ) lcnt++;
					
					for ( int j = 0 ; j < links[1][c].m_num ; ++j )
					{
						int t = links[1][c][j];
						if ( level[t] == N )
						{
							if ( mark.notexist(t) )
							{
								mark.insert(t,2);
								q.push_back(t);
							}
							else
							{
								if ( mark.get(t) == 1 )
								{
									printf("DAG Property Error!!!\n");
								}
							}
						}
					}
				}
			}
			q.clean();
		
			level[p] = id;			

			l2n[id] = p;

			mark.clean();
		}

		tstop = clock();

		printf( "Time: %0.2lf\n",(double)(tstop-tstart)/CLOCKS_PER_SEC );

		printf( "Label Size: %0.0lf\n",lcnt*4.0);
		
		pn.free_mem();
	}

/*	void load_hl_label( string st )
	{
		mark.clean();

		string labelfn = "hl\\"+st+".HRTindexO";

		ifstream in(labelfn);

		string buf, inbuf, sub;

		int vid;
        getline(in,buf);
                for(int i = 0; i < N; i++) {
                        getline(in,buf);
                        buf = buf.substr(buf.find_first_of(":")+2);
                        // parse Lin_l
                        inbuf = buf.substr(0,buf.find_first_of("#"));
                        while (inbuf.find(" ")!=string::npos) {
                                sub = inbuf.substr(0,inbuf.find(" "));
                                istringstream is(sub);
								is >> vid;
								labels[1][i].push_back(vid);
								if ( mark.notexist(vid) )
								{
									mark.insert(vid,1);									
								}
								else
								{
									mark.inc(vid,1);
								}
                                inbuf.erase(0,inbuf.find(" ")+1);
                        }
                        // parse Lout_l
                        buf.erase(0, buf.find_first_of("#")+2);
                        while (buf.find(" ")!=string::npos) {
                                sub = buf.substr(0,buf.find(" "));
                                istringstream(sub) >> vid;
								labels[0][i].push_back(vid);
								if ( mark.notexist(vid) )
								{
									mark.insert(vid,1);									
								}
								else
								{
									mark.inc(vid,1);
								}
                                buf.erase(0,buf.find(" ")+1);
								if ( buf[0] == '#' ) break;
                        }
                }


		q.clean();

		for ( int i = 0 ; i < N ; ++i )
		{
			if ( mark.get(i) == 2 )
			{
				q.push_back(i);
			}
		}

		for ( int i = 0 ; i < q.m_num ; ++i )
		{
			int p = q[i];
			level[p] = N-1-i;

			for ( dir = 0 ; dir < 2 ; ++dir )
			{
				for ( int j = 0 ; j < labels[dir][p].m_num ; ++j )
				{
					int c = labels[dir][p][j];
					mark.dec( c );
					if ( mark.get( c ) == 2 ) q.push_back(c);
				}
			}
		}

		if ( q.m_num != N )
		{
			int t = N-q.m_num;
			for ( int i  = 0 ; i < N ; ++i )
			{
				if ( mark.get(i) > 2 )
				{
					t--;
					level[i] = t;
				}
			}
		}

		mark.clean();
		q.clean();
	}
*/
	void load_tf_level( string st )
	{
		int x, y = 0;

		//load level
		string levelfn = "index\\"+st+"_level";
		FILE *Lfile = fopen(levelfn.c_str(),"r");

		mark.clean();

		fscanf(Lfile,"%d",&x);

		for ( int i = 0 ; i < x ; ++i )
		{
			fscanf(Lfile,"%d",&y);
			level[y] = i;
			l2n[i] = y;
			mark.insert(y,-2);
		}

		for ( int i = 0 ; i < N ; ++i )
		{
			if ( fscanf(Lfile,"%d",&y) != 1 ) printf("Error!!!\n");
			if ( mark.notexist(i) ) mark.insert(i,y);
		}

		fclose(Lfile);

		y = N;
		
		for ( int r = 0 ; r < N ; ++r )
		{
			bool flag = true;

			for ( int i = N-1 ; i >= 0 ; --i )
			{
				if ( mark.get(i) == r )
				{
					--y;
					flag = false;
					level[i] = y;
					l2n[y] = i;
				}
			}

			if ( flag )
			{
				for ( int i = N-1 ; i >= 0 ; --i )
				{
					if ( mark.get(i) == -1 )
					{
						--y;
						level[i] = y;
						l2n[y] = i;
					}
				}
				break;
			}
		}

		if ( x != y ) printf("Error!!!\n");

		mark.clean();

		iVector<Key_Value<double,int> > pn;

		//compute hierarchies
		lcnt = 0;

		for ( int i = 0 ; i < N ; ++i )
		{
			pn.push_back( Key_Value<double,int>( -1.0*(links[0][i].m_num+1)*(links[1][i].m_num+1), i ) );
			level[i] = N;
		}

		pn.Sort();

		for ( int i = TOPNUM-1 ; i >= 0 ; --i )
		{
			level[pn[i].value] = i-TOPNUM;
		}
	}

	void compute_index_tf( int* rlevel )
	{
		tstart = clock();

		iVector<Key_Value<int,int> > pn;

		//compute hierarchies
		lcnt = 0;

		for ( int i = 0 ; i < N ; ++i )
		{
			pn.push_back( Key_Value<int,int>( rlevel[i], i ) );
			level[i] = N;
		}

		pn.Sort();

		for ( int id = 0 ; id < N ; ++id )
		{
			int p = pn[id].value;			

			//find children of p, mark them, and create back labels
			q.push_back(p);
			mark.insert(p,2);

			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];
				for ( int j = 0 ; j < links[0][c].m_num ; ++j )
				{
					int t = links[0][c][j];
					if ( level[t] == N )
					{
						if ( mark.notexist(t) )
						{
							mark.insert(t,1);
							lcnt += AddLabel( t, p, 1, id );							
							q.push_back(t);		
						}
					}
				}
			}
			AddLabel( p, p, 0, id );
			q.clean();

			//find ancestors of p
			q.push_back(p);
			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];
				for ( int j = 0 ; j < links[1][c].m_num ; ++j )
				{
					int t = links[1][c][j];
					if ( level[t] == N )
					{
						if ( mark.notexist(t) )
						{
							mark.insert(t,2);
							lcnt += AddLabel( t, p, 0, id );
							q.push_back(t);
						}
						else
						{
							if ( mark.get(t) == 1 )
							{
								printf("DAG Property Error!!!\n");
							}
						}
					}
				}
			}
			AddLabel( p, p, 1, id );
			q.clean();
		
			level[p] = id;			

			l2n[id] = p;

			mark.clean();
		}

		tstop = clock();

		printf( "Time: %0.2lf\n",(double)(tstop-tstart)/CLOCKS_PER_SEC );

		printf( "Label Size: %0.0lf\n",lcnt*4.0);
		
		pn.free_mem();
	}

	void free()
	{
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < N ; ++i )
			{
				labels[dir][i].free_mem();
				links[dir][i].free_mem();
			}
			labels[dir].free_mem();
			links[dir].free_mem();
		}

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			cost[dir].free_mem();
		}
		pq.m_data.free_mem();
		pq.pos.free_mem();
		q.free_mem();
		mark.free_mem();
		//random_nodes.free_mem();
		delete[] sum;
	}

	void generate_subgraph( string graph, int rsize, iVector<int>* neighbors, iVector<int>* nn )
	{
		random_nodes.initialize(N);

		for ( int i = 0 ; i < rsize ; )
		{
			int x = rand()<<15;
			x += rand();
			x %= N;
		//	int x = N-rsize+i;
			if ( random_nodes.notexist(x) )
			{
				random_nodes.insert( x, N-rsize+i );				
				nodes.push_back( x );
				i++;
			}
		}

		for ( int i = 0, j = 0 ; i < N ; ++i )
		{
			if ( random_nodes.notexist(i) )
			{
				random_nodes.insert( i , j );
				++j;
			}
		}

		if( nodes.m_num != rsize ) printf("Random nodes generating error!!!");

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < rsize ; ++i )
			{
				int p = nodes[i];
				int xx = 0;
				for ( int j = 0; j < links[dir][p].m_num ; ++j )
				{
					int c = links[dir][p][j];
					if ( random_nodes.get(c) < N-rsize+i ) 
					{
						neighbors[dir].push_back( random_nodes.get(c) );
						xx++;
					}
				}
				nn[dir].push_back(xx);
			}
		}

		string subgraph = "SB"+graph;
		string oq = "NT"+graph;
		//string oq = "1.txt";
		string sq = "NTSB"+graph;

		FILE* file = fopen(subgraph.c_str(),"w");

		fprintf(file,"%d\n",N-rsize);

		for ( int i = 0 ; i < N ; ++i )
		{
			if ( random_nodes.get(i)+rsize < N )
			{
				level[random_nodes.get(i)] = level[i];
				fprintf(file,"%d:",random_nodes.get(i));
				for ( int j = 0 ; j < links[0][i].m_num ; ++j )
				{
					int t = links[0][i][j];
					if ( random_nodes.get(t)+rsize < N )
						fprintf(file," %d",random_nodes.get(t));
				}
				fprintf(file," -1\n");
			}
		}

		fclose(file);

		FILE *qfile = fopen( oq.c_str(), "r" );
		FILE *sqfile = fopen( sq.c_str(), "w" );
		
		int x,y;
		for ( ; fscanf(qfile,"%d %d",&x,&y)==2 ; )
		{
			fprintf(sqfile,"%d %d\n",random_nodes.get(x),random_nodes.get(y));
		}

		fclose(qfile);
		fclose(sqfile);
	}

	void generate_DG( string graph, int rsize )
	{
		random_nodes.initialize(N);

		for ( int i = 0 ; i < rsize ; )
		{
			int x = rand()<<15;
			x += rand();
			x %= N;
		//	int x = N-rsize+i;
			if ( random_nodes.notexist(x) )
			{
				random_nodes.insert( x, N-rsize+i );				
				nodes.push_back( x );
				i++;
			}
		}

		for ( int i = 0, j = 0 ; i < N ; ++i )
		{
			if ( random_nodes.notexist(i) )
			{
				random_nodes.insert( i , j );
				++j;
			}
		}

		if( nodes.m_num != rsize ) printf("Random nodes generating error!!!");

		string ups = "UP"+graph;

		FILE *ofile = fopen(ups.c_str(),"w");

		for ( int i = 0 ; i < rsize ; ++i )
		{
			int p = nodes[i];

			fprintf( ofile, "a %d ",N-rsize+i );

			for ( dir = 0 ; dir < 2 ; ++dir )
			{
				int cnt = 0;

				for ( int j = 0 ; j < links[1-dir][p].m_num ; ++j )
				{
					int c = links[1-dir][p][j];
					if ( random_nodes.get(c) < N-rsize+i )
					{
						cnt++;
					}
				}

				fprintf( ofile, "%d ", cnt );

				for ( int j = 0 ; j < links[1-dir][p].m_num ; ++j )
				{
					int c = links[1-dir][p][j];
					if ( random_nodes.get(c) < N-rsize+i )
					{
						fprintf(ofile,"%d ",random_nodes.get(c));
					}
				}
			}

			fprintf( ofile,"\n");
		}

		for ( int i = 0 ; i < rsize ; ++i )
		{
			fprintf(ofile, "r %d\n", N-rsize+i);
		}

		fclose(ofile);

		string subgraph = "P"+graph;

		FILE* file = fopen(subgraph.c_str(),"w");

		fprintf(file,"graph_for_greach\n%d\n",N-rsize);

		for ( int i = 0 ; i < N ; ++i )
		{
			if ( random_nodes.get(i)+rsize < N )
			{
				level[random_nodes.get(i)] = level[i];
				fprintf(file,"%d: ",random_nodes.get(i));
				for ( int j = 0 ; j < links[0][i].m_num ; ++j )
				{
					int t = links[0][i][j];
					if ( random_nodes.get(t)+rsize < N )
						fprintf(file,"%d ",random_nodes.get(t));
				}
				fprintf(file,"#\n");
			}
		}

		fclose(file);
	}

	void generate_Grail( string graph )
	{
		string Grail_graph = "G_"+graph;

		FILE* file = fopen(Grail_graph.c_str(),"w");

		fprintf(file,"graph_for_greach\n%d\n",N);

		for ( int i = 0 ; i < N ; ++i )
		{
			fprintf(file,"%d: ",i);
			for ( int j = 0 ; j < links[0][i].m_num ; ++j )
			{
				int t = links[0][i][j];
				fprintf(file,"%d ",t);				
			}
			fprintf(file,"#\n");
		}

		fclose(file);
	}

	void AddNode( int x, iVector<int>* neighbors, int* ss, iVector<int>* nn, int px )
	{
		//compute level

		//if ( x == 128553 )
		//{
		//	printf("fuck\n");
		//}

		//compute potential candis
		candi.clean();
		mark.clean();
		list[0].clean(); list[1].clean();

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = ss[dir] ; i < ss[dir]+nn[dir][px] ; ++i )
			{
				int p = neighbors[dir][i];
				for ( int j = 0 ; j < labels[dir][p].m_num ; ++j )
				{
					int c = labels[dir][p][j];
					if ( mark.notexist( c ) ) 
					{
						if ( levelhash.nottop(TOPNUM,c) ) 
						{
							candi.push_back( Triple<int>( levelhash.N2L[c], c, dir ) );
							mark.insert( c, c );
						}
						else
						{
							lcnt += AddLabel1(x,c,dir);
							mark.insert( c, -10 );
							list[dir].push_back( c );
						}
					}
				}
			}			
		}

		//compute level

		//if ( x == 598648 )
		//{
		//	printf("fuck\n");
		//}

		//clean candis while generating left and right wings
		if ( candi.m_num > 0 ) candi.Sort();
		
		int j = 0;

		for ( int i = 0 ; i < candi.m_num ; ++i )
		{
			int p = candi[i].y;

			if ( mark.get(p) == p )
			{
				if ( i != j )
				{
					candi[j] = candi[i];
				}
				j++;

				int d = candi[i].w;

				mark.insert( p, 0 );

				for ( int k = 0 ; k < backlabels[1-d][p].m_num ; ++k )
				{
					int c = backlabels[1-d][p][k];

					if ( c != p && ( mark.notexist(c) || mark.get(c) == c ) )
					{
						mark.insert( c, p );
						mark.inc( p );
					}
				}
			}
		}

		candi.m_num = j;

		//scan candis from lowest level to highest level
		int sum = 0, k = 0, y = levelhash.Last;

		n2x[0].clean(); n2x[1].clean();

		for ( int i = candi.m_num-1 ; i >= 0 ; --i )
		{
			int c = candi[i].y;
			int l = candi[i].x;
			int d = candi[i].w;

			//if ( x == 128553 )
			//{
			//	printf("c: %d, d: %d, son: %d\n",c,d,mark.get(c));
			//}

			for ( int j = 0 ; j < backlabels[1-d][c].m_num ; ++j )
			{
				int t = backlabels[1-d][c][j];

				if ( mark.get(t) == c && t != c && !FuckCheck2( t, x, 1-d, levelhash ) )
				{
					sum++;
				}
			}
			for ( int j = 0 ; j < n2x[1-d].m_num ; ++j )
			{
				int t = n2x[1-d][j];
				if ( !FuckCheck( t, c, d, x ) )
				{
					sum--;
				}
			}

			if ( sum < k )
			{
				k = sum;
				y = l;
			}			

			for ( int j = 0 ; j < backlabels[1-d][c].m_num ; ++j )
			{
				int t = backlabels[1-d][c][j];
				if ( t != c && mark.get(t) == c )
				{
					n2x[d].push_back(t);
				}
			}
			n2x[d].push_back(c);

			mark.insert( c, d-2 );
		}

		//label update
		for ( int i = 0 ; i < mark.occur.m_num ; ++i )
		{
			int p = mark.occur[i];

			if ( mark.get(p) >= 0 )
			{
				int c = mark.get(p);

				if ( mark.get(c) >= 0 )
				{
					printf("Error in marking for non crutial nodes\n");
				}
				candi.push_back(Triple<int>(levelhash.N2L[c],p,mark.get(c)+2));
			}
		}

		if ( candi.m_num > 0 ) candi.Sort();

		for ( int i = 0 ; i < candi.m_num ; ++i )
		{
			int p = candi[i].y;
			int d = candi[i].w;

			if ( mark.get(p) < 0 )
			{
				if ( levelhash.N2L[p] < y )
				{
					lcnt += AddLabel1(x,p,d);
				}
				else
				{
					lcnt += AddLabel1(p,x,1-d);
					for ( int k = 0 ; k < backlabels[d][x].m_num ; ++k )
					{
						if ( backlabels[d][x][k] != p ) cleanmark.insert( backlabels[d][x][k], 1 );
					}
				}
				list[d].push_back( p );
			}
			else if ( candi[i].x >= y )
			{
					lcnt += AddLabel1(p,x,1-d);
					for ( int k = 0 ; k < backlabels[d][x].m_num ; ++k )
					{
						if ( backlabels[d][x][k] != p ) cleanmark.insert( backlabels[d][x][k], 1 );
					}
			}
			
			for ( int j = 0 ; j < list[1-d].m_num ; ++j )
			{
				int c = list[1-d][j];

				if ( AddLabel1(p,c,1-d) == 1 )
				{
					lcnt++;
					for ( int k = 0 ; k < backlabels[d][c].m_num ; ++k )
					{
						if ( backlabels[d][c][k] != c && backlabels[d][c][k] != p ) cleanmark.insert( backlabels[d][c][k], 1 );
					}
				}
			}

			if ( cleanmark.occur.m_num > 0 )
			{
				int k = 0;
				for ( int j = 0 ; j < labels[1-d][p].m_num ; ++j )
				{
					if ( cleanmark.notexist( labels[1-d][p][j] ) )
					{
						if ( k != j )
						{
							labels[1-d][p][k] = labels[1-d][p][j];
						}
						k++;
					}
					else 
					{
						if ( !backlabels[1-d][labels[1-d][p][j]].remove(p) )
						{
							printf("Backlabel deleting Error!!!\n");
						}
						lcnt--;
					}
				}
				labels[1-d][p].m_num = k;

				cleanmark.clean();
			}
		}

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			lcnt += AddLabel1(x,x,dir);
		}

		//level update
		levelhash.insert( x, y );
	}

	void compute_index_p()
	{
		tstart = clock();

//		father.initialize(N);

		//compute topological order and cost
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < N ; ++i )
			{
				sum[i] = links[1-dir][i].m_num;
				if ( sum[i] == 0 ) q.push_back(i);
				cost[dir][i] = 1;
			}		

			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int p = q[i];

				for ( int j = 0 ; j < links[dir][p].m_num ; ++j )
				{
					int t = links[dir][p][j];
					cost[dir][t] += cost[dir][p]/links[dir][p].m_num;
					sum[t]--;
					if ( sum[t] == 0 ) q.push_back(t);
				}
			}

			q.clean();
		}

		//compute hierarchies
		lcnt = 0;

		double *tcost = new double[N];

		for ( int i = 0 ; i < N ; ++i )
		{
			pq.insert(i,-(cost[0][i]-1)*(cost[1][i]-1)/(cost[0][i]+cost[1][i]));
			level[i] = N;
		}

		for ( int id = 0 ; id < N ; ++id )
		{
			//get top node with lazy update
			int p ;

			for ( ; true ; )
			{
				p = pq.head();
				pq.insert( p, -(cost[0][p]-1)*(cost[1][p]-1)/(cost[0][p]+cost[1][p]) );
				if ( pq.head() == p )
				{
					p = pq.pop();
					break;
				}
			}

		//	father.clean();

			//find children of p, mark them, and create back labels
			q.push_back(p);
			mark.insert(p,2);
			tcost[p] = 0;

		//	if ( id % 10000 == 0 ) printf("%d AVG LAB NUM: %0.2lf\n",id,lcnt*1.0/N);

			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];
				if ( AddLabel( c, p, 1, id ) )
				{
					if ( c != p ) lcnt++;
					for ( int j = 0 ; j < links[0][c].m_num ; ++j )
					{
						int t = links[0][c][j];
						if ( level[t] == N )
						{
							if ( mark.notexist(t) )
							{
								mark.insert(t,1);
								q.push_back(t);		
								tcost[t] = 0;
							}
						}
					}
				}
			}
			q.clean();

			//find ancestors of p
			q.push_back(p);
			for ( int i = 0 ; i < q.m_num ; ++i )
			{
				int c = q[i];

				if ( AddLabel( c, p, 0, id ) )
				{
					if ( c != p ) lcnt++;
					for ( int j = 0 ; j < links[1][c].m_num ; ++j )
					{
						int t = links[1][c][j];
						if ( level[t] == N )
						{
							if ( mark.notexist(t) )
							{
								mark.insert(t,2);
								q.push_back(t);
								tcost[t] = 0;
							}
							else
							{
								if ( mark.get(t) == 1 )
								{
									printf("DAG Property Error!!!\n");
								}
							}
						}
					}
				}
			}
			q.clean();

			for ( int i = 0 ; i < mark.occur.m_num ; ++i )
			{
				int c = mark.occur[i];
				if ( mark.get(c) == 2 )
				{
					for ( int j = 0 ; j < links[0][c].m_num ; ++j )
					{
						int t = links[0][c][j];
						if ( level[t] == N && ( t == p || mark.get(t)==1 ) )
						{
							tcost[c] += cost[1][t]/links[1][t].m_num;
							tcost[t] += cost[0][c]/links[0][c].m_num;
						}
					}
				}
			}

			for ( int i = 1 ; i < mark.occur.m_num ; ++i )
			{
				int c = mark.occur[i];
				int d = mark.get(c)-1;
				cost[d][c] -= tcost[c];

				for ( int j = 0 ; j < links[d][c].m_num ; ++j )
				{
					int t = links[d][c][j];

					if ( level[t] == N )
					{
						tcost[t] += tcost[c]/links[d][c].m_num;
					}
				}
			}

			level[p] = id;			

			l2n[id] = p;

			mark.clean();
		}

		reduce();

		tstop = clock();

		printf( "Time: %0.2lf\n",(double)(tstop-tstart)/CLOCKS_PER_SEC );

		printf( "Label Size: %0.0lf\n",lcnt*4.0);

		delete[] tcost;
	}



	void AppendNode( int x, iVector<int>* neighbors, int* ss, iVector<int>* nn, int px )
	{
		//compute level

		//if ( x == 128553 )
		//{
		//	printf("fuck\n");
		//}

		//compute potential candis
		candi.clean();
		mark.clean();
		list[0].clean(); list[1].clean();

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = ss[dir] ; i < ss[dir]+nn[dir][px] ; ++i )
			{
				int p = neighbors[dir][i];
				for ( int j = 0 ; j < labels[dir][p].m_num ; ++j )
				{
					int c = labels[dir][p][j];
					if ( mark.notexist( c ) ) 
					{
						if ( levelhash.nottop(TOPNUM,c) ) 
						{
							candi.push_back( Triple<int>( levelhash.N2L[c], c, dir ) );
							mark.insert( c, c );
						}
						else
						{
							lcnt += AddLabel1(x,c,dir);
							mark.insert( c, -10 );
							list[dir].push_back( c );
						}
					}
				}
			}			
		}

		//compute level

		//if ( x == 598648 )
		//{
		//	printf("fuck\n");
		//}

		//clean candis while generating left and right wings
		if ( candi.m_num > 0 ) candi.Sort();
		
		int j = 0;

		for ( int i = 0 ; i < candi.m_num ; ++i )
		{
			int p = candi[i].y;

			if ( mark.get(p) == p )
			{
				if ( i != j )
				{
					candi[j] = candi[i];
				}
				j++;

				int d = candi[i].w;

				mark.insert( p, 0 );

				for ( int k = 0 ; k < backlabels[1-d][p].m_num ; ++k )
				{
					int c = backlabels[1-d][p][k];

					if ( c != p && ( mark.notexist(c) || mark.get(c) == c ) )
					{
						mark.insert( c, p );
						mark.inc( p );
					}
				}
			}
		}

		candi.m_num = j;

		int y = levelhash.Last;

		for ( int i = candi.m_num-1 ; i >= 0 ; --i )
		{
			int c = candi[i].y;
			int l = candi[i].x;
			int d = candi[i].w;

			mark.insert( c, d-2 );
		}

		//label update
		for ( int i = 0 ; i < mark.occur.m_num ; ++i )
		{
			int p = mark.occur[i];

			if ( mark.get(p) >= 0 )
			{
				int c = mark.get(p);

				if ( mark.get(c) >= 0 )
				{
					printf("Error in marking for non crutial nodes\n");
				}
				candi.push_back(Triple<int>(levelhash.N2L[c],p,mark.get(c)+2));
			}
		}

		if ( candi.m_num > 0 ) candi.Sort();

		for ( int i = 0 ; i < candi.m_num ; ++i )
		{
			int p = candi[i].y;
			int d = candi[i].w;

			if ( mark.get(p) < 0 )
			{
				if ( levelhash.N2L[p] < y )
				{
					lcnt += AddLabel1(x,p,d);
				}
				else
				{
					lcnt += AddLabel1(p,x,1-d);
					for ( int k = 0 ; k < backlabels[d][x].m_num ; ++k )
					{
						if ( backlabels[d][x][k] != p ) cleanmark.insert( backlabels[d][x][k], 1 );
					}
				}
				list[d].push_back( p );
			}
			else if ( candi[i].x >= y )
			{
					lcnt += AddLabel1(p,x,1-d);
					for ( int k = 0 ; k < backlabels[d][x].m_num ; ++k )
					{
						if ( backlabels[d][x][k] != p ) cleanmark.insert( backlabels[d][x][k], 1 );
					}
			}
			
			for ( int j = 0 ; j < list[1-d].m_num ; ++j )
			{
				int c = list[1-d][j];

				if ( AddLabel1(p,c,1-d) == 1 )
				{
					lcnt++;
					for ( int k = 0 ; k < backlabels[d][c].m_num ; ++k )
					{
						if ( backlabels[d][c][k] != c && backlabels[d][c][k] != p ) cleanmark.insert( backlabels[d][c][k], 1 );
					}
				}
			}

			if ( cleanmark.occur.m_num > 0 )
			{
				int k = 0;
				for ( int j = 0 ; j < labels[1-d][p].m_num ; ++j )
				{
					if ( cleanmark.notexist( labels[1-d][p][j] ) )
					{
						if ( k != j )
						{
							labels[1-d][p][k] = labels[1-d][p][j];
						}
						k++;
					}
					else 
					{
						if ( !backlabels[1-d][labels[1-d][p][j]].remove(p) )
						{
							printf("Backlabel deleting Error!!!\n");
						}
						lcnt--;
					}
				}
				labels[1-d][p].m_num = k;

				cleanmark.clean();
			}
		}

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			lcnt += AddLabel1(x,x,dir);
		}

		//level update
		levelhash.insert( x, y );
	}

	void AddDelTest( int rsize, vector<int>* tests)
	{
		srand((unsigned int)time(NULL));
		levelhash.initialize( level.m_data, N, 10 );

		modified.re_allocate(N); modified.m_num = N;

		opq.initialize( N );

		random_nodes.initialize(N);

		for ( int i = 0 ; i < rsize ; )
		{
			int x;
			if(tests->size()>0)  x = (*tests)[i];
			else  x = rand()%N;//<<15;
//			x += rand();
//			x %= N;
			if ( random_nodes.notexist(x) && levelhash.nottop(TOPNUM,x) )
			{
				random_nodes.insert( x, N-rsize+i );				
				nodes.push_back( x );
				i++;
			}
		}

		//for ( int i = 0 ; i < links[1][110769].m_num ; ++i )
		//{
		//	if ( random_nodes.exist(links[1][110769][i]) )
		//	{
		//		printf("inlink: %d\n",links[1][110769][i]);
		//	}
		//}

		if( nodes.m_num != rsize ) printf("Random nodes generating error!!!");

		iVector<int> neighbors[2],nn[2];

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < rsize ; ++i )
			{
				int p = nodes[i];
				int xx = 0;
				for ( int j = 0; j < links[dir][p].m_num ; ++j )
				{
					int c = links[dir][p][j];
					if ( random_nodes.notexist(c) || random_nodes.get(c) < N-rsize+i ) 
					{
						neighbors[dir].push_back( c );
						xx++;
					}
				}
				nn[dir].push_back(xx);
			}
		}

		printf("\n");

		clock_t tstart,tstop;

		tstart = clock();

		for ( int i = 0 ; i < rsize ; ++i )
		{
			// cout<< "deleting: " << i << endl;
			DeleteNode( nodes[i] );
			// if ( (i+1) % 100 == 0 )
			// {
				// printf("\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r%d nodes deleted",i+1);
			// }
		}
		printf("\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\rDeletion Done. ");
		tstop = clock();

		printf( " Time: %0.2lf ms\n",(double)(tstop-tstart)*1000/CLOCKS_PER_SEC );

		printf("Label Size: %0.0lf\n",lcnt*4.0);

		int s[2]={0,0};

		tstart = clock();

		for ( int i = 0 ; i < rsize ; ++i )
		{
			// cout<< "inserting: " << i << endl;
			AddNode( nodes[i], neighbors, s, nn, i );
			s[0] += nn[0][i];
			s[1] += nn[1][i];		

			// if ( (i+1) % 100 == 0 )
			// {
			// 	printf("\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r%d nodes inserted",i+1);
			// }
		}
		printf("\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\r\rInsertion Done. ");
		tstop = clock();

		printf( " Time: %0.5lf ms\n",(double)(tstop-tstart) * 1000.0 /CLOCKS_PER_SEC );

		printf("Label Size: %0.0lf\n",lcnt*4.0);
	}

	bool FuckCheck( int s , int t , int d, int z )
	{
		int i,j,x,y;

		for ( i = 0, j = 0 ; i < labels[d][s].m_num && j < labels[1-d][t].m_num ; )
		{
			x = labels[d][s][i]; y = labels[1-d][t][j];

			if ( x == y && x != t && x != z ) return true;
			else if ( x < y ) ++i;
			else ++j;
		}

		return false;
	}

	bool FuckCheck2( int s , int t , int d, LevelHash& lh )
	{
		int i,j,x,y;

		for ( i = 0, j = 0 ; i < labels[d][s].m_num && j < labels[1-d][t].m_num ; )
		{
			x = labels[d][s][i]; y = labels[1-d][t][j];

			if ( x == y && !lh.nottop(TOPNUM,x) ) return true;
			else if ( x < y ) ++i;
			else ++j;
		}

		return false;
	}

	int UpgradeNode( int x )
	{
		int r = lcnt;
		//compute level

		//if ( x == 2921696 )
		//{
		//	printf("fuck\n");
		//}

		//compute potential candis
		candi.clean();
		mark.clean();
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < links[dir][x].m_num ; ++i )
			{
				int p = links[dir][x][i];
				for ( int j = 0 ; j < labels[dir][p].m_num ; ++j )
				{
					int c = labels[dir][p][j];
					if ( mark.notexist( c ) && levelhash.nottop(TOPNUM, c) )
					{
						mark.insert( c, c );
						candi.push_back( Triple<int>( levelhash.N2L[c], c, dir ) );
					}
				}
			}			
		}

		//clean candis while generating left and right wings
		if ( candi.m_num > 0 ) candi.Sort();
		
		int j = 0;

		for ( int i = 0 ; i < candi.m_num ; ++i )
		{
			int p = candi[i].y;

			if ( mark.get(p) == p )
			{
				if ( i != j )
				{
					candi[j] = candi[i];
				}
				j++;

				int d = candi[i].w;

				mark.insert( p, 0 );

				for ( int k = 0 ; k < backlabels[1-d][p].m_num ; ++k )
				{
					int c = backlabels[1-d][p][k];

					if ( c != p && ( mark.notexist(c) || mark.get(c) == c ) )
					{
						mark.insert( c, p );
						mark.inc( p );
					}
				}
			}
		}

		candi.m_num = j;

		//scan candis from lowest level to highest level
		int sum = 0, k = 0, y = levelhash.Last, tmp = 0;

		n2x[0].clean(); n2x[1].clean();

		bool flag = true;

		for ( int i = candi.m_num-1 ; i >= 0 ; --i )
		{
			int c = candi[i].y;
			int l = candi[i].x;
			int d = candi[i].w;

			if ( l < levelhash.N2L[x] && flag )
			{
				flag = false;
				tmp = sum;
			}

			//if ( x == 2921696 )
			//{
			//	if ( flag ) printf("low: ");
			//	printf("c: %d, d: %d, son: %d\n",c,d,mark.get(c));
			//}

//			sum += mark.get(c);
			for ( int j = 0 ; j < backlabels[1-d][c].m_num ; ++j )
			{
				int t = backlabels[1-d][c][j];

				if ( mark.get(t) == c && t != c && !FuckCheck2( t, x, 1-d, levelhash ) )
				{
					sum++;
				}
			}
			for ( int j = 0 ; j < n2x[1-d].m_num ; ++j )
			{
				int t = n2x[1-d][j];
				if ( !FuckCheck( t, c, d, x ) )
				{
					sum--;
				}
			}

			if ( sum < k )
			{
				k = sum;
				y = l;
			}			

			for ( int j = 0 ; j < backlabels[1-d][c].m_num ; ++j )
			{
				int t = backlabels[1-d][c][j];
				if ( t != c && mark.get(t) == c )
				{
					n2x[d].push_back(t);
				}
			}
			n2x[d].push_back(c);

			mark.insert( c, d-2 );
		}

		if ( flag ) tmp = sum;
		
		if ( tmp == k ) return 0;

		OPR++;

		//label update
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			int j = 0;
			for ( int i = 0 ; i < labels[dir][x].m_num ; ++i )
			{
				int c = labels[dir][x][i];
				if ( levelhash.nottop(TOPNUM,c) )
				{
					if ( c != x )
					{
						if ( !backlabels[dir][labels[dir][x][i]].remove( x ) )
						{
							printf("Backlabel deleting Error!!!\n");
						}
						lcnt--;
					}
				}
				else
				{
					labels[dir][x][j] = labels[dir][x][i];
					j++;
				}
			}
			labels[dir][x].m_num = j;

			for ( int i = 0 ; i < backlabels[dir][x].m_num ; ++i )
			{
				if ( backlabels[dir][x][i] != x )
				{					
					if ( !labels[dir][backlabels[dir][x][i]].remove(x) )
					{
						printf("Label deleting Error!!!\n");
					}
					lcnt--;
				}
			}
			backlabels[dir][x].clean();
		}

		for ( int i = 0 ; i < mark.occur.m_num ; ++i )
		{
			int p = mark.occur[i];

			if ( mark.get(p) >= 0 )
			{
				int c = mark.get(p);

				if ( mark.get(c) >= 0 )
				{
					printf("Error in marking for non crutial nodes\n");
				}
				candi.push_back(Triple<int>(levelhash.N2L[c],p,mark.get(c)+2));
			}
		}

		if ( candi.m_num > 0 ) candi.Sort();
		list[0].clean(); list[1].clean();

		for ( int i = 0 ; i < candi.m_num ; ++i )
		{
			int p = candi[i].y;
			int d = candi[i].w;

			if ( mark.get(p) < 0 )
			{
				if ( levelhash.N2L[p] < y )
				{
					lcnt += AddLabel1(x,p,d);
				}
				else
				{
					lcnt += AddLabel1(p,x,1-d);
					for ( int k = 0 ; k < backlabels[d][x].m_num ; ++k )
					{
						if ( backlabels[d][x][k] != p ) cleanmark.insert( backlabels[d][x][k], 1 );
					}
				}
				list[d].push_back( p );
			}
			else if ( candi[i].x >= y )
			{
					lcnt += AddLabel1(p,x,1-d);
					for ( int k = 0 ; k < backlabels[d][x].m_num ; ++k )
					{
						if ( backlabels[d][x][k] != p ) cleanmark.insert( backlabels[d][x][k], 1 );
					}
			}
			
			for ( int j = 0 ; j < list[1-d].m_num ; ++j )
			{
				int c = list[1-d][j];

				if ( AddLabel1(p,c,1-d) == 1 )
				{
					lcnt++;
					for ( int k = 0 ; k < backlabels[d][c].m_num ; ++k )
					{
						if ( backlabels[d][c][k] != c && backlabels[d][c][k] != p ) cleanmark.insert( backlabels[d][c][k], 1 );
					}
				}
			}

			if ( cleanmark.occur.m_num > 0 )
			{
				int k = 0;
				for ( int j = 0 ; j < labels[1-d][p].m_num ; ++j )
				{
					if ( cleanmark.notexist( labels[1-d][p][j] ) )
					{
						if ( k != j )
						{
							labels[1-d][p][k] = labels[1-d][p][j];
						}
						k++;
					}
					else 
					{
						if ( !backlabels[1-d][labels[1-d][p][j]].remove(p) )
						{
							printf("Backlabel deleting Error!!!\n");
						}
						lcnt--;
					}
				}
				labels[1-d][p].m_num = k;

				cleanmark.clean();
			}
		}

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			AddLabel1(x,x,dir);
		}

		//level update
		levelhash.swap( x, y );


		if ( r-lcnt != tmp-k )
		{
			if ( candi.m_num <= 100 )
			{
				printf("fuck %d %d\n",x,candi.m_num);
			}
		}
		return tmp-k;
	}

//       void build_backlink( int extra )
//        {
//                iVector<int> tmp;
//                for ( dir = 0 ; dir < 2 ; ++dir )
//                {
//                        mark.clean();
//                        for ( int i = 0 ; i < N ; ++i )
//                        {
//                                tmp.clean();
//                                tmp.push_back( labels[dir][i].m_data, labels[dir][i].m_num );
//                                labels[dir][i].clean();
//                                labels[dir][i].re_allocate( tmp.m_num+20 );
//                                labels[dir][i].push_back( tmp.m_data, tmp.m_num );
//                                for ( int j = 0 ; j < labels[dir][i].m_num ; ++j )
//                                {
//                                        int c = labels[dir][i][j];
//                                        if ( mark.notexist(c) ) mark.insert( c, 1 );
//                                        else mark.inc(c);
//                                }
//                        }
////              }
//
////              for ( dir = 0 ; dir < 2 ; ++dir )
////              {
//                        backlabels[dir].re_allocate( N+extra );
//                        backlabels[dir].m_num = N;
//
//                        for ( int i = 0 ; i < N ; ++i )
//                        {
//                                backlabels[dir][i].re_allocate( mark.get(i)+20 );
//                        }
//
//
//                        for ( int i = 0 ; i < N ; ++i )
//                        {
//                                for ( int j = 0 ; j < labels[dir][i].m_num ; ++j )
//                                {
//                                        labels[dir][i][j] =     l2n[labels[dir][i][j]];
//                                        backlabels[dir][labels[dir][i][j]].push_back( i );
//                                }
//                                if ( labels[dir][i].m_num > 0 ) labels[dir][i].Sort();
//                        }
//                }
//        }


	void build_backlink( int extra )
	{
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			backlabels[dir].re_allocate( N+extra );
			backlabels[dir].m_num = N;
			for ( int i = 0 ; i < N ; ++i )
			{
				for ( int j = 0 ; j < labels[dir][i].m_num ; ++j )
				{
					labels[dir][i][j] =	l2n[labels[dir][i][j]];
					backlabels[dir][labels[dir][i][j]].push_back( i );
				}			
				if ( labels[dir][i].m_num > 0 ) labels[dir][i].Sort();
			}
		}
	}


	void collapse( int s, int t, int *color )
	{
		q.clean();
		mark.clean();

		q.push_back(s);

		for ( int i = 0 ; i < q.m_num ; ++i )
		{
			int p = q[i];

			if ( p != s )
			{
				for ( dir = 0 ; dir < 2 ; ++dir )
				{
					labels[dir][s].push_back( labels[dir][p].m_data, labels[dir][p].m_num );
					labels[dir][p].clean();
				}
				color[p] = color[s];
			}

			for ( int j = 0 ; j < links[0][p].m_num ; ++j )
			{
				int c = links[0][p][j];

				if ( reachable( c, t ) )
				{
					if ( mark.notexist(c) )
					{
						mark.insert(c,p);
						q.push_back(c);
					}
				}
			}
		}

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			labels[dir][s].unique();
			for ( int i = 0 ; i < labels[dir][s].m_num ; ++i )
			{
				int p = labels[dir][s][i],k=0;
				int flag = -1;

				for ( int j = 0 ; j < backlabels[dir][p].m_num ; ++j )
				{
					int c = backlabels[dir][p][j];
					if ( color[c] == color[s] )
					{
						if ( flag > -1 ) continue;
						else
						{
							flag = k;
							c = s;
						}
					}
					backlabels[dir][p][k] = c;					
					k++;
				}

				backlabels[dir][p].m_num = k;

				for ( ; flag < k-1 && flag > 0 ; )
				{
					if ( backlabels[dir][p][flag] > backlabels[dir][p][flag+1] )
					{
						int tmp = backlabels[dir][p][flag];
						backlabels[dir][p][flag] = backlabels[dir][p][flag+1];
						backlabels[dir][p][flag+1] = tmp;
					}
					else if ( backlabels[dir][p][flag] < backlabels[dir][p][flag-1] )
					{
						int tmp = backlabels[dir][p][flag];
						backlabels[dir][p][flag] = backlabels[dir][p][flag-1];
						backlabels[dir][p][flag-1] = tmp;
					}
					else break;
				}
			}
		}		
	}

	bool TriangleReduceCheck( int s , int t , int d )
	{
		int i,j,x,y;

		for ( i = 0, j = 0 ; i < labels[d][s].m_num && j < labels[1-d][t].m_num ; )
		{
			x = labels[d][s][i]; y = labels[1-d][t][j];

			if ( x == y && x != t && x != s ) 
			{
				return true;
			}
			else if ( x < y ) ++i;
			else ++j;
		}

		return false;
	}

	bool check( int s, int t, int d )
	{
		int i,j,x,y;

		for ( i = 0, j = 0 ; i < backlabels[1-d][s].m_num && j < backlabels[d][t].m_num ; )
		{
			x = backlabels[1-d][s][i]; y = backlabels[d][t][j];

			if ( x == y )
			{
				if ( x != s && x != t ) return true;
				i++;
				j++;
			}
			else if ( x < y ) ++i;
			else ++j;
		}

		return false;
	}

	void DeleteNode( int x )
	{
		//update label
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			opq.clean();
			opq.insert( x, order[1-dir][x] );
			mark.clean();
			mark.insert(x,-1);

			for ( ; !opq.empty() ; )
			{
				int p = opq.pop();

				modified[p] = false;

				bool flag = false;

				if ( p != x )
				{
					flag = true;
					for ( int k = 0 ; k < links[dir][p].m_num ; ++k )
					{
						if ( modified[links[dir][p][k]] )
						{
							flag = false;
							break;
						}
					}
				}
				
				if ( !flag )
				{
					rnodes.clean();

					if ( p != x )
					{
						//mark reachable nodes
						//rnodes.clean();
						for ( int k = 0 ; k < links[dir][p].m_num ; ++k )
						{
							int t = links[dir][p][k];

							if ( t != x && t != p )
							{
								for ( int h = 0 ; h < labels[dir][t].m_num ; ++h )
								{
									int c = labels[dir][t][h];

									if ( mark.get(c) != p )
									{
										mark.insert( c , p );
										rnodes.push_back( c );
									}
								}
							}
						}
					}
				
					//remove unreachable nodes from labels[dir][p]
					int j = 0;
					for ( int k = 0 ; k < labels[dir][p].m_num ; ++k )
					{
						int t = labels[dir][p][k];

						if ( ( t != p || t == x ) && mark.get(t) != p )
						{
							backlabels[dir][t].remove(p);
							lcnt--;
							modified[p] = true;
						}
						else
						{
							mark.insert(t,-1);
							labels[dir][p][j] = t;
							j++;
						}
					}
					labels[dir][p].m_num = j;
					//add newly reachable nodes to labels[dir][p]
					for ( int k = 0 ; k < rnodes.m_num ; ++k )
					{
						int t = rnodes[k];

						if ( mark.get(t) == p && levelhash.N2L[t] < levelhash.N2L[p] )
						{
							if ( AddLabel1(p,t,dir) ) 
							{
								modified[p] = true;
								lcnt++;
							}
						}
					}
				}

				//extend p
				if ( levelhash.nottop(TOPNUM,p) )
				{
					for ( int j = 0 ; j < links[1-dir][p].m_num ; ++j )
					{
						int c = links[1-dir][p][j];

						if ( mark.notexist(c) )
						{
							opq.insert(c,order[1-dir][c]);
							mark.insert(c, -1);
						}
					}
				}
			}
		}

		//links update
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < links[dir][x].m_num ; ++i )
			{
				int c = links[dir][x][i];

				links[1-dir][c].remove_unsorted(x);
			}
			links[dir][x].clean();
		}


		//level update
		levelhash.remove(x);
	}

	void testDel()
	{
		levelhash.initialize( level.m_data, N, 1 );

		modified.re_allocate(N); modified.m_num = N;

		opq.initialize( N );

		random_nodes.initialize( N );

		srand(time(0));

		for ( int i = 0 ; i < 10000 ; )
		{
		//	if ( i%100 == 0 ) printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d ",i);

			//if ( nodes[i] == 20912 )
			//{
			//	printf("fuck\n");
			//}

			int x = rand();
			x = (x<<15)+rand();
			x %= N;

			if ( random_nodes.notexist( x ) )
			{
				random_nodes.insert( x , 1 );
				DeleteNode( x );
				i++;
			}
		}

		printf("Label size after deletion: %d\n",lcnt*4);
	}

	bool reachable( int s, int t )
	{
		int i , j , x , y;
		for ( i = 0 , j = 0 ; i < labels[0][s].m_num && j < labels[1][t].m_num ; )
		{
			x = labels[0][s][i]; y = labels[1][t][j];
			if ( x == y ) 
			{
				return true;
			}
			else if ( x < y ) i++;
			else j++;
		}
		return false;
	}

	void outputlinks( int x, int d )
	{
		printf("%d:",x);
		for ( int i = 0 ; i < links[d][x].m_num ; ++i )
		{
			int c = links[d][x][i];

			printf(" %d",c);
		}
		printf("\n");
	}

	void compute_order()
	{
		int *sum = new int[N];
		for ( dir = 0 ; dir < 2 ; ++dir )//dir: in-edges and out edges
		{
			q.clean();
			order[dir].re_allocate(N);
			order[dir].m_num = N;

			for ( int i = 0 ; i < N ; ++i )
			{
				sum[i] = links[1-dir][i].m_num;//leaf nodes
				if ( sum[i] == 0 ) q.push_back(i);
				order[dir][i] = 0;
			}		

			for ( int i = 0 ; i < q.m_num ; ++i ) //start at leaf nodes?
			{
				int p = q[i];

				for ( int j = 0 ; j < links[dir][p].m_num ; ++j )
				{
					int t = links[dir][p][j];
					order[dir][t] = max( order[dir][t], order[dir][p]+1 );
					sum[t]--;
					if ( sum[t] == 0 ) q.push_back(t);
				}
			}

			q.clean();
		}

		//iVector<Key_Value<int,int> > on;

		//for ( dir = 0 ; dir < 2 ; ++dir )
		//{
		//	for ( int i = 0 ; i < N ; ++i )
		//	{
		//		on.clean();
		//		for ( int j = 0 ; j < links[dir][i].m_num ; ++j )
		//		{
		//			int c = links[dir][i][j];
		//			on.push_back(Key_Value<int,int>( order[dir][c], c ) );
		//		}
		//		if ( on.m_num > 0 ) on.Sort();
		//		for ( int j = 0 ; j < on.m_num ; ++j )
		//		{
		//			links[dir][i][j] = on[j].value;
		//		}
		//	}
		//}
	}

	void AddEdge( int* V )
	{
		if ( reachable( V[0] , V[1] ) ) return;

		list[1].clean(); list[0].clean();
		
		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < labels[1-dir][V[dir]].m_num ; ++i )
			{
				int p = labels[1-dir][V[dir]][i];
				list[1-dir].push_back( p );

				for ( int j = 0 ; j < backlabels[1-dir][V[1-dir]].m_num ; ++i )
				{
					lcnt += AddLabel1( backlabels[1-dir][V[1-dir]][j], p, 1-dir );
				}		
			}
		}

		for ( int i = 0 ; i < list[1].m_num ; ++i )
		{
			int s = list[1][i];

			for ( int j = 0 ; j < list[0].m_num ; ++j )
			{
				int t = list[0][j];

				if ( level[s] < level[t] )
				{
					for ( int k = 0 ; k < backlabels[1][t].m_num ; ++k )
					{
						int tt = backlabels[1][t][k];
						lcnt += AddLabel1( tt, s, 1 );
					}
				}
				else
				{
					for ( int k = 0 ; k < backlabels[0][s].m_num ; ++k )
					{
						int ss = backlabels[0][s][k];
						lcnt += AddLabel1( ss, t, 0 );
					}
				}
			}
		}
	}


	void reduce()
	{
		//for ( dir = 0 ; dir < 2 ; ++dir )
		//{
		//	backlabels[dir].re_allocate(N);
		//	backlabels[dir].m_num = N;
		//	for ( int i = 0 ; i < N ; ++i )
		//	{
		//		int p = l2n[i];

		//		for ( int j = 0 ; j < labels[dir][p].m_num ; ++j )
		//		{
		//			backlabels[dir][labels[dir][p][j]].push_back(i);					
		//		}
		//	}
		//}

		int m = 0;

		for ( dir = 0 ; dir < 2 ; ++dir )
		{
			for ( int i = 0 ; i < N ; ++i )
			{
				int x = 0;
				for ( int j = 0 ; j < labels[dir][i].m_num ; ++j )
				{
					int p = labels[dir][i][j];

					int q = l2n[p];					

					bool flag = true;

					if ( q != i )
					{
						for ( int k = 0 ; k < labels[1-dir][q].m_num ; ++k )
						{
							if ( mark.get( labels[1-dir][q][k] ) == i )
							{
								lcnt--;
								flag = false;
								break;
							}
						}
					}

					if ( flag )
					{
						labels[dir][i][x] = p;
						x++;
					}

					mark.insert( p, i );
				}
				labels[dir][i].m_num = x;
			}
		}

		printf("%0.0lf\n",lcnt*4.0);
	}
 

	void generateRG( FILE* file, int n, int D )
	{
		links[0].re_allocate(n);

		for ( int i = 0 ; i < n*D ; ++i )
		{
			int x = (rand()<<15);
			x += rand();
			x %= n;

			int y = (rand()<<15);
			y += rand();
			y %= n;

			if ( x < y )
			{
				links[0][x].push_back(y);
			}
			else if ( y < x )
			{
				links[0][y].push_back(x);
			}
		}

		for ( int i = 0 ; i < n ; ++i )
		{
			links[0][i].unique();
		}

		fprintf(file, "%d\n", n );

		double sum = 0;

		for ( int i = 0 ; i < n ; ++i )
		{
			if ( links[0][i].m_num > 0 )
			{
				fprintf(file, "%d: ", i);

				for ( int j = 0 ; j < links[0][i].m_num ; ++j )
				{
					fprintf(file, "%d ",links[0][i][j]);
					sum++;
				}

				fprintf(file, "-1\n");
			}
		}

		printf("Real Degree: %0.2lf\n",sum/n);
	}
};



int main(int argc, char** argv)
{
	string GraphFileName, IndexFileStem, Updates;

	GraphFileName = ""; IndexFileStem = ""; Updates="";

	int style=0;
	int opr = 0;

	for ( int i = 1 ; i < argc ; ++i )
	{
		if ( strcmp( argv[i], "-g" ) == 0 ) 
		{
			GraphFileName = string(argv[i+1]);
			++i;
		}
		
		if ( strcmp( argv[i], "-upd" ) == 0 )
		{
			opr = 0;
			Updates = string(argv[i+1]);
			// cout<< Updates << endl;
		}
	}

	Butterfly run;

	run.load_graph( GraphFileName, 0);

	if ( style == 0 )
	{
		run.compute_order();
		run.compute_index_r1();
	}
	else if ( style == 1 )
	{
		run.compute_index_p1();
	}
	run.build_backlink(0);
	if ( opr == 0 )
	{
		run.compute_order();
		vector<int> tests;
		std::ifstream file(Updates);
		if (!file.is_open()) {
			std::cerr << "random testing" << std::endl;
			run.AddDelTest(10000, &tests);
		}
		else{
			string line;
			while (std::getline(file, line)) {
				if (!line.empty()) {
					int number = stoi(line);
					tests.push_back(number);
				}
			}
			file.close(); // 关闭文件
			cout << "Updating ... \n" << endl;
			run.AddDelTest(tests.size(), &tests);
		}
	}

	return 0;
}

