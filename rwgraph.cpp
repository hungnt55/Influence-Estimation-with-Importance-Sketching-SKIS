#include "rwgraph.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <cstdlib>
#include <unistd.h>
#include <sstream>

using namespace std;

double Graph::getTotalProb()
{
	return totalProb;
}

const vector<int> & Graph::operator [] (int u) const
{
	return adjList[u];
}


const vector<int> & Graph::operator [] (int u)
{
	return adjList[u];
}


const vector<UI> & Graph::getWeight (int u) const
{
        return weights[u];
}

const vector<UI> & Graph::getWeight (int u)
{
        return weights[u];
}

/*
* get degree of node u
*/
int Graph::getDegree(int u) const
{
	return adjList[u].size();
}

/*
* get the number of nodes
*/
int Graph::getSize() const
{
	return numNodes;
}

/*
* get the number of edges
*/
int Graph::getEdge() const
{
	return numEdges;
}

void Graph::findEdge(int edgeId, int & end1, int & end2)
{
//	int curNode = 1;
//	while (node_deg[curNode] < 2 || adjIndex[curNode][node_deg[curNode]-2] < edgeId){
//		curNode++;
//	}

//	end1 = curNode;
//	end2 = adjList[curNode][edgeId - adjIndex[curNode][0]];
}

/*
* read binary graph input for LT model
* difference between LT and IC is for LT we accumulate the weights for fast choosing a random node
*/
void Graph::readGraphLT(const char* filename, float scale)
{
   	FILE * pFile;
    	pFile = fopen(filename, "rb");
    	fread(&numNodes, sizeof(int), 1, pFile);
    	fread(&numEdges, sizeof(long long), 1, pFile);
    	node_deg=vector<int>(numNodes + 1);
    	fread(&node_deg[1], sizeof(int), numNodes, pFile);
        
	vector<int> a;
    	vector<UI> b;
    	adjList.push_back(a);
    	weights.push_back(b);
//	adjIndex.push_back(a);
	int counter = 0;
        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<int> tmp(node_deg[i]);
		vector<int> tmp2(node_deg[i]);
                fread(&tmp[0], sizeof(int), node_deg[i], pFile);

                adjList.push_back(tmp);
		for (int j = 0; j < node_deg[i]; ++j){
			tmp2[j] = counter;
			counter++;
		}
//		adjIndex.push_back(tmp2);

        }

        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<float> tmp(node_deg[i] + 1, 0);
                vector<UI> tmp1(node_deg[i] + 1, 0);
                fread(&tmp[1], sizeof(float), node_deg[i], pFile);

                for(int j = 1;j < node_deg[i] + 1; ++j){
                        tmp[j] += tmp[j-1];
			tmp1[j] = (UI)scale*tmp[j]*UI_MAX;
		}
		
		tmp1[node_deg[i]] = UI_MAX;
		
                weights.push_back(tmp1);
		node_deg[i]++;
        }
}

/*
* read input graph for IC model
*/
void Graph::readGraphIC(const char* filename)
{
    	FILE * pFile;
    	pFile = fopen(filename, "rb");
    	fread(&numNodes, sizeof(int), 1, pFile);
    	fread(&numEdges, sizeof(long long), 1, pFile);
    	node_deg=vector<int>(numNodes + 1);
    	fread(&node_deg[1], sizeof(int), numNodes, pFile);
        
	vector<int> a;
	vector<UI> b;
    	adjList.push_back(a);
    	weights.push_back(b);
	
        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<int> tmp(node_deg[i]);
                fread(&tmp[0], sizeof(int), node_deg[i], pFile);
                adjList.push_back(tmp);
        }

	totalProb = 0;
	nodeProb = vector<UI>(numNodes+1,0);
	vector<double> nodeProbT(numNodes+1,0);
        for (unsigned int i = 1; i <= numNodes; ++i){
                vector<float> tmp(node_deg[i] + 1, 0);
		vector<UI> tmp1(node_deg[i] + 1, 0);
                fread(&tmp[1], sizeof(float), node_deg[i], pFile);

		double multiplier = 1;
                for(int j = 1;j < node_deg[i] + 1; ++j){
                        tmp1[j] = (tmp[j]*multiplier)*UI_MAX + tmp1[j-1];
			multiplier *= (1-tmp[j]);
                }
		totalProb += (1-multiplier);
		nodeProbT[i] = totalProb;

		if (tmp1[node_deg[i]] <= 0)
                        tmp1[node_deg[i]] = UI_MAX;
		
                weights.push_back(tmp1);
        }
	for (unsigned int i = 1; i <= numNodes; ++i){
		nodeProb[i] = (nodeProbT[i]/totalProb)*UI_MAX;
	}
	
	vector<double>().swap(nodeProbT);
}

void Graph::writeToFile(const char * filename)
{/*
	ofstream output(filename);
	for (unsigned int i = 0; i < numNodes; ++i){
		for (unsigned int j = 0; j < adjList[i].size(); ++j){
			if (adjList[i][j] > i){
				output << adjList[i][j] << " " << i << " " << weights[i][j] << endl;
			}
		}
	}
	output.close();
*/	
}

// choose a random edge in LT model based on linear search
inline int HyperGraph::randIndex_lin(const vector<UI> &w, unsigned int si, sfmt_t &sfmtSeed)
{
        UI ranNum = sfmt_genrand_uint32(&sfmtSeed);
        if (si <= 1 || ranNum > w[si - 1])
                return -1;

        for (unsigned int i = 1; i < si; ++i){
                if (ranNum <= w[i])
                        return i;
        }
        return -1;
}

// choose a random live edge in LT model based on binary search
/*inline int HyperGraph::randIndex_bin(const vector<UI> &w, unsigned int si, sfmt_t &sfmtSeed)
{
	UI ran = sfmt_genrand_uint32(&sfmtSeed);
	if (si <= 1 || ran > w[si - 1])
                return -1;
        int left = 1;
        int right = si - 1;
        int prob;
        for (unsigned int i = 0; i < si; ++i){
                prob = (left + right)/2;
                if (w[prob - 1] > ran){
                        right = prob - 1;
                        continue;
                }
                if (w[prob] <= ran){
                        left = prob + 1;
                        continue;
                }
                break;
        }
        return prob;
}*/

inline int HyperGraph::randIndex_bin(const vector<UI> &w, unsigned int si, sfmt_t &sfmtSeed, int ben, int end, bool normalized)
{
        UI ran = sfmt_genrand_uint32(&sfmtSeed);
        if (normalized){
                ran = ran*(w[end-1]/(double)UI_MAX);
        } else {
                ran = ran*((UI_MAX-w[ben])/(double)UI_MAX) + w[ben];
        }

        if (si <= 1 || ran > w[end - 1] || ben >= end-1)
                return -1;
        int left = ben;
        int right = end - 1;
        int prob = ben;
        for (unsigned int i = 0; i < si; ++i){
                prob = (left + right)/2;
		if (w[prob] <= ran){
                        left = prob + 1;
                        continue;
                }
                if (prob-1 >= ben && w[prob - 1] > ran){
                        right = prob - 1;
                        continue;
                }
                break;
        }
        return prob;
}

HyperGraph::HyperGraph(unsigned int n, unsigned int m)
{
	sfmt_init_gen_rand(&sfmtSeed, rand());
	node_edge = vector<vector<int> >(n+1);
	node_root = vector<int>(n+1,0);
	maxDegree = 0;
	numNodes = n;
	numEdges = m;
	curEdge=0;
}

void HyperGraph::updateDeg(){
	unsigned int num=edge_node.size();
	for (unsigned int i = curEdge; i < num; ++i){
		unsigned int num2 = edge_node[i].size();
		for (unsigned int j=0;j<num2;++j){
			node_edge[edge_node[i][j]].push_back(i);
		}
	}
	curEdge = edge_node.size();
	edgeMask = vector<bool> (num, false);
	visitEdge = vector<int> (num, 0);
	num_masked = 0;
}

void HyperGraph::updateEdge(){
	curEdge = edge_node.size();
}

/*
* Add a hyperedge into the hypergraph
*/
void HyperGraph::addEdge(vector<int> & edge)
{
	edge_node.push_back(edge);
	unsigned int ind = edge_node.size() - 1;
	for (unsigned int i = 0; i < edge.size(); ++i)
		node_edge[edge[i]].push_back(ind);
}

/*
* Add a hyperedge into the hypergraph while keeping track of the node with max degree
*/
void HyperGraph::addEdgeD(vector<int> & edge)
{
        edge_node.push_back(edge);
        int ind = edge_node.size() - 1;
        for (unsigned int i = 0; i < edge.size(); ++i){
                node_edge[edge[i]].push_back(ind);
		if (node_edge[edge[i]].size() > maxDegree)
			maxDegree = node_edge[edge[i]].size();
	}
}

/*
* get an edge from the hypergraph
*/
const vector<int> & HyperGraph::getEdge(int e) const{
	return edge_node[e];
}

const vector<int> & HyperGraph::getEdge(int e){
	return edge_node[e];
}

/*
* get the list of hyperedges incident to node n
*/
const vector<int> & HyperGraph::getNode(int n) const{
	return node_edge[n];
}

const vector<int> & HyperGraph::getNode(int n){
	return node_edge[n];
}

/*
* get the number of hyperedges
*/
int HyperGraph::getNumEdge() const
{
        return edge_node.size();
}

/*
* get the maximum degree
*/
int HyperGraph::getMaxDegree()
{
	return maxDegree;
}

/*
* remove all the hyperedges
*/
void HyperGraph::clearEdges()
{
	edge_node.clear();
	node_edge.clear();
	cout << "clear edges!" << endl;
       maxDegree = 0;
}

void HyperGraph::setSeeds(vector<int> &s)
{
	seeds = vector<bool>(numNodes+1, false);
	for (unsigned int i = 0; i < s.size(); ++i){
		seeds[s[i]] = true;
	}
	
	for (unsigned int i = 1; i < numNodes; ++i){
		activeNodes.push_back(i);
	}

	numNodes = activeNodes.size();
}

/*
* polling process under LT model
*/ 
int HyperGraph::pollingLT2(Graph &g, vector<unsigned int> & link, vector<bool> & interdict, bool & success, unsigned int k, vector<bool> &visit, vector<int> &visit_mark, vector<double> & prob, unsigned int & num_marked, sfmt_t & sfmtSeed1)
{	
	unsigned int i;
        unsigned int gSize = g.getSize();
        unsigned int cur = activeNodes[sfmt_genrand_uint32(&sfmtSeed1)%numNodes];
        num_marked = 0;
	bool hit = false;
        for (i = 0; i < gSize; ++i){
                if (visit[cur] == true) break;

		if (interdict[cur]){
                        success = true;
                }

                visit[cur] = true;
                visit_mark[num_marked] = cur;
		num_marked++;
		if (link[cur]){
                        if (rand() <= prob[cur]*RAND_MAX){
                                hit = true;
                                break;
                        }
                }

                int ind;
                //if (g.weights[cur].size() >= 32)
//                        ind = randIndex_bin(g.weights[cur],g.node_deg[cur],sfmtSeed1);
                //else
                        ind = randIndex_lin(g.weights[cur],g.node_deg[cur],sfmtSeed1);

                if (ind == -1)
                        break;

                cur = g.adjList[cur][ind - 1];
        }

	if (hit){
		return 1;
	}
	success = false;
	return 0;
}

int HyperGraph::pollingLT3(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark, vector<double> & prob, unsigned int & num_marked, sfmt_t & sfmtSeed1)
{
        unsigned int i;
        unsigned int gSize = g.getSize();
        unsigned int cur = activeNodes[sfmt_genrand_uint32(&sfmtSeed1)%numNodes];
        num_marked = 0;
        bool hit = false;
        for (i = 0; i < gSize; ++i){
                if (visit[cur] == true) break;

                visit[cur] = true;
                visit_mark[num_marked] = cur;
                num_marked++;
                if (link[cur]){
                        if (rand() <= prob[cur]*RAND_MAX){
                                hit = true;
                                break;
                        }
                }

                int ind;
                //if (g.weights[cur].size() >= 32)
    //                    ind = randIndex_bin(g.weights[cur],g.node_deg[cur],sfmtSeed1);
                //else
                        ind = randIndex_lin(g.weights[cur],g.node_deg[cur],sfmtSeed1);

                if (ind == -1)
                        break;

                cur = g.adjList[cur][ind - 1];
        }

        if (hit){
                return 1;
        }

        return 0;
}

void HyperGraph::pushHyperedges(vector<vector<int> > & hyperedges){
	for (unsigned int i = 0; i < hyperedges.size(); ++i){
		edge_node.push_back(hyperedges[i]);
	}
}

bool HyperGraph::pollingLT(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{
        unsigned int i;
        bool t = false;
        unsigned int gSize = g.getSize();
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%gSize;
        unsigned int num_marked = 0;
        for (i = 0; i < gSize; ++i){
		if (link[cur] < k){
                        t=true;
			break;
                }
                if (visit[cur] == true) break;
                visit[cur] = true;
                visit_mark[num_marked] = cur;
		num_marked++;
		int ind;
		//if (g.weights[cur].size() >= 32)
  //                      ind = randIndex_bin(g.weights[cur],g.node_deg[cur],sfmtSeed);
                //else
                        ind = randIndex_lin(g.weights[cur],g.node_deg[cur],sfmtSeed);

                if (ind == -1)
                        break;
                cur = g.adjList[cur][ind - 1];
        }
        for (i = 0; i < num_marked; ++i){
                visit[visit_mark[i]]=false;
        }
        return t;
}


void HyperGraph::pollingLT1(Graph &g, vector<bool> &visit, vector<int> &visit_mark, vector<bool> & link, vector<double> & prob)
{
	unsigned int i;
        unsigned int gSize = g.getSize();
        unsigned int cur = activeNodes[sfmt_genrand_uint32(&sfmtSeed)%numNodes];
        unsigned int num_marked = 0;
        bool hit = false;
        for (i = 0; i < gSize; ++i){
                if (visit[cur] == true) break;
                visit[cur] = true;
                visit_mark[num_marked] = cur;
                num_marked++;
		if (link[cur]){
                        if (rand() <= prob[cur]*RAND_MAX){
                                hit = true;
                                break;
                        }
                }

                int ind;
                //if (g.weights[cur].size() >= 32)
    //                    ind = randIndex_bin(g.weights[cur],g.node_deg[cur],sfmtSeed);
                //else
                        ind = randIndex_lin(g.weights[cur],g.node_deg[cur],sfmtSeed);

                if (ind == -1)
                        break;
                cur = g.adjList[cur][ind - 1];
        }
        if (hit){
                edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
//		cout << num_marked << endl;
        }

        for (i = 0; i < num_marked; ++i){
                visit[visit_mark[i]]=false;
        }
}

void HyperGraph::pollingIC1(Graph &g, vector<bool> &visit, vector<int> &visit_mark)
{
        int i;
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%(g.getSize());
        int num_marked=1;
	int curPos=0;
        visit[cur] = true;
        visit_mark[0] = cur;
        while(curPos < num_marked){
                cur = visit_mark[curPos];
		curPos++;
                const vector<UI> &w=g.getWeight(cur);
                const vector<int> &neigh = g[cur];
                for (i = 0; i < g.node_deg[cur]; ++i){
                        if (sfmt_genrand_uint32(&sfmtSeed) <  w[i+1]){
                        	if (!visit[neigh[i]]){
                                        visit[neigh[i]] = true;
                                        visit_mark[num_marked]=neigh[i];
					num_marked++;
                                }
                        }
                }
        }
	edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
        for(i = 0; i < num_marked;++i){
                visit[visit_mark[i]]=false;
        }
}

/*
* polling process under IC model

int HyperGraph::pollingIC2(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark)
{
        int i;
        unsigned int cur = sfmt_genrand_uint32(&sfmtSeed)%(g.getSize());
	int curPos=0;
        int num_marked=1;
	visit[cur] = true;
	visit_mark[0] = cur;
	unsigned int t = k;
        while(curPos < num_marked){
		cur = visit_mark[curPos];
		curPos++;
		if (link[cur] < t)
	               t=link[cur];
		const vector<UI> &w=g.getWeight(cur);
                const vector<int> &neigh = g[cur];
		for (i = 0; i < g.node_deg[cur]; ++i){
			if (sfmt_genrand_uint32(&sfmtSeed) <  w[i+1]){
				if (!visit[neigh[i]]){
					visit[neigh[i]] = true;
					visit_mark[num_marked]=neigh[i];
					num_marked++;
				}
			}
		}
        }
        edge_node.push_back(vector<int>(visit_mark.begin()+1,visit_mark.begin()+num_marked));

	for(i = 0; i < num_marked;++i){
		visit[visit_mark[i]]=false;
	}
	return t;
}
*/

void HyperGraph::pollingIC(Graph &g, vector<bool> &visit, vector<int> &visit_mark, unsigned int & num_marked, sfmt_t & sfmtSeed1)
{
        int selectedPos;
	unsigned int cur = randIndex_bin(g.nodeProb,g.numNodes+1,sfmtSeed1,0,g.numNodes+1,false);
        num_marked = 0;
       	visit[cur] = true;
        visit_mark[0] = cur;
       	num_marked++;

	const vector<UI> &w=g.getWeight(cur);
       	int wsize = w.size();
        const vector<int> &neigh = g[cur];
       	selectedPos = randIndex_bin(w,wsize,sfmtSeed1,0,wsize,true);
        while (selectedPos > -1){
       	        if (!visit[neigh[selectedPos-1]]){
               	        visit[neigh[selectedPos-1]] = true;
                       	visit_mark[num_marked]=neigh[selectedPos-1];
                        num_marked++;
       	         }
               	 selectedPos = randIndex_bin(w,wsize,sfmtSeed1,selectedPos,wsize,false);
       	}

        unsigned int curPos=1;
       	while(curPos < num_marked){
                cur = visit_mark[curPos];
       	        curPos++;
               	const vector<UI> &w=g.getWeight(cur);
                int wsize = w.size();
       	        const vector<int> &neigh = g[cur];
               	selectedPos = randIndex_bin(w,wsize,sfmtSeed1,0,wsize,false);
                while (selectedPos > -1){
       	                if (!visit[neigh[selectedPos-1]]){
               	                visit[neigh[selectedPos-1]] = true;
                       	        visit_mark[num_marked]=neigh[selectedPos-1];
                               	num_marked++;
                        }
       	                selectedPos = randIndex_bin(w,wsize,sfmtSeed1,selectedPos,wsize,false);
               	}
       	}

       	for(unsigned int i = 0; i < num_marked;++i){
               	visit[visit_mark[i]]=false;
       	}
	#pragma omp critical
	{
		edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
		node_root[visit_mark[0]]++;
	}
}

int HyperGraph::estimateInf(vector<int> & seeds){
	num_masked = 0;
	for (int i : seeds){
		for (int j : node_edge[i]){
			if (!edgeMask[j]){
				edgeMask[j] = true;
				visitEdge[num_masked] = j;
				num_masked++;
			}
		}
	}
	
	for (int i = 0; i < num_masked; ++i){
		edgeMask[visitEdge[i]] = false;
	}
	for (int i : seeds){
		num_masked -= node_root[i];
	}
	return num_masked;
}

int HyperGraph::pollingIC2(Graph &g, vector<unsigned int> & link, unsigned int k, vector<bool> &visit, vector<int> &visit_mark, unsigned int & num_marked, sfmt_t & sfmtSeed1)
{
        int selectedPos;
        unsigned int cur = randIndex_bin(g.nodeProb,g.numNodes+1,sfmtSeed1,0,g.numNodes+1,true);
        num_marked=0;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;
	unsigned int t = g.numNodes+1;

        const vector<UI> &w=g.getWeight(cur);
        int wsize = w.size();
        const vector<int> &neigh = g[cur];
        selectedPos = randIndex_bin(w,wsize,sfmtSeed1,0,wsize,true);
        while (selectedPos > -1){
                if (!visit[neigh[selectedPos-1]]){
                        visit[neigh[selectedPos-1]] = true;
                        visit_mark[num_marked]=neigh[selectedPos-1];
                        num_marked++;
                 }
                 selectedPos = randIndex_bin(w,wsize,sfmtSeed1,selectedPos,wsize,false);
        }

        unsigned int curPos=1;
        while(curPos < num_marked){
                cur = visit_mark[curPos];
		if (link[cur] > 0)
                       t=link[cur];
                curPos++;
                const vector<UI> &w=g.getWeight(cur);
                int wsize = w.size();
                const vector<int> &neigh = g[cur];
                selectedPos = randIndex_bin(w,wsize,sfmtSeed1,0,wsize,false);
                while (selectedPos > -1){
                        if (!visit[neigh[selectedPos-1]]){
                                visit[neigh[selectedPos-1]] = true;
                                visit_mark[num_marked]=neigh[selectedPos-1];
                                num_marked++;
                        }
                        selectedPos = randIndex_bin(w,wsize,sfmtSeed1,selectedPos,wsize,false);
                }
        }

        for(unsigned int i = 0; i < num_marked;++i){
                visit[visit_mark[i]]=false;
        }
	//edge_node.push_back(vector<int>(visit_mark.begin(),visit_mark.begin()+num_marked));
	if (t <= g.numNodes)
		return 1;
	else 
		return 0;
}

/*
* convert from an integer to a string
*/
string intToStr(int i) {
        stringstream ss;
        ss << i;
        return ss.str();
}

/*
* convert from a strong to an integer
*/
unsigned int strToInt(string s) {
        unsigned int i;
        istringstream myStream(s);

        if (myStream>>i) {
                return i;
        } else {
                cout << "String " << s << " is not a number." << endl;
                return atoi(s.c_str());
        }
        return i;
}

/*
* measure the consumed memory
*/
float getCurrentMemoryUsage() {

        string pid = intToStr(unsigned(getpid()));
        string outfile = "tmp_" + pid + ".txt";
        string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;
        system(command.c_str());

        string mem_str;
        ifstream ifs(outfile.c_str());
        std::getline(ifs, mem_str);
        ifs.close();

        mem_str = mem_str.substr(0, mem_str.size()-1);
        float mem = (float)strToInt(mem_str);

	command = "rm " + outfile;
        system(command.c_str());

        return mem/1024;

        return 0;
}
