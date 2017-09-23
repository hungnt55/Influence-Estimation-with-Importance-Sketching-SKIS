#ifndef _HYPERGRAPH_H_
#define _HYPERGRAPH_H_

#include "rwgraph.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <random>
#include "mappedHeap.hpp"
#include "HeapData.hpp"

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_set_num_threads(int t) { return 1;}
inline omp_int_t omp_get_thread_num() { return 0;}
#endif

using namespace std;

/*
* building the hypergraph procedure which generates hyperedges following LT model
*/
long long addHyperedge(Graph & g, HyperGraph & hg, int t, long long num, bool lt, vector<int> & seeds, vector<double> & prob)
{
	int numNodes = g.getSize();

//	std::uniform_int_distribution<int> intGen(1,numNodes);

	// temporary array for holding hyperedges before put them into the hypergraph
//	vector<vector<vector<int> > > tmp(t);

	omp_set_num_threads(t);
	
	long long iter = 0;
	int c = 10000;

        #pragma omp parallel
	{
	//	int id = omp_get_thread_num();
//                int ver;
	
		vector<int> visit_mark(numNodes+1,0);
		vector<bool> visit(numNodes+1,false);
		vector<bool> link(numNodes+1,false);
		for (unsigned int i = 0; i < seeds.size(); ++i){
			link[seeds[i]] = true;
		}

		while (iter < num){
       	                for (int i = 0; i < c; ++i){
			//	vector<int> he;
	                       	hg.pollingLT1(g,visit,visit_mark,link,prob);

//				hg.addEdge(he);
                               	//tmp[id].push_back(he);
                       	}
		
			#pragma omp atomic
			iter += c;
		}
	}
	//long long numE =0;
	// add hyperedges to hypergraph
	/*
	for (int i = 0; i < t; ++i){
		numE += tmp[i].size();
		for (long long j = 0; j < tmp[i].size(); ++j)
                       hg.addEdge(tmp[i][j]);
		// clean memory after taking all the edges
		vector<vector<int> >().swap(tmp[i]);
	}*/
	hg.updateDeg();
	return hg.getNumEdge();
}

/*
* find seed nodes procedure using greedy algorithm
*/
int buildSeedSet(HyperGraph & hg, vector<int> & seeds, unsigned int n, int k, vector<double> &degree, double & ob)
{	
	long long i;
	unsigned int j,l,maxInd;
	vector<int> e, nList;

	vector<int> nodeDegree(n,0);
	vector<int> indx(n,0);
	for (j = 0; j < n; ++j){
		indx[j] = j;
		nodeDegree[j] = hg.getNode(j).size();
	}

	InfCost<int> hd(&nodeDegree[0]);
	MappedHeap<InfCost<int> > heap(indx,hd);
	long long numEdge = hg.getNumEdge();
	vector<bool> edgeMark(numEdge, false);
	vector<bool> nodeMark(n+1, true);
	double totalCost = 0;
	i=1;
	int deg = 0;
	int leftover = 0;	

	while(totalCost < k && !heap.empty()){
		maxInd = heap.pop();
		nodeMark[maxInd] = false;
		totalCost++;
		e = hg.getNode(maxInd);
		degree[i] = degree[i-1]+nodeDegree[maxInd];
		deg += nodeDegree[maxInd];
		seeds.push_back(maxInd);
		for (j = 0; j < e.size(); ++j){
			if (edgeMark[e[j]]){
				continue;
			}
			nList = hg.getEdge(e[j]);
			for (l = 0; l < nList.size(); ++l){
				nodeDegree[nList[l]]--;
				if (nodeMark[nList[l]]){
					heap.heapify(nList[l]);
				}
			}
			edgeMark[e[j]] = true;
		}
		i++;
	}
	
	for (i = 0; i < k; ++i){
		maxInd = heap.pop();
                leftover += nodeDegree[maxInd];
	}
	vector<int>().swap(nodeDegree);
	vector<int>().swap(e);
	vector<int>().swap(nList);
	vector<int>().swap(indx);
	vector<bool>().swap(edgeMark);
	ob = deg/(deg+leftover);
	return deg;
}

int buildSeedSetLinear(HyperGraph & hg, vector<int> & seeds, unsigned int n, int k, vector<double> &degree, double & ob)
{
        long long i;
        unsigned int j,l,maxInd, maxDeg;
        vector<int> e, nList;

        vector<unsigned int> nodeDegree(n,0);
        for (j = 0; j < n; ++j){
                nodeDegree[j] = hg.getNode(j).size();
	//	cout << j << " " << nodeDegree[j] << endl;
        }

        long long numEdge = hg.getNumEdge();
        vector<bool> edgeMark(numEdge, false);
        vector<bool> nodeMark(n+1, true);
        double totalCost = 0;
        i=1;
	int deg = 0;
	ob = 0;
        while(totalCost < k){
                maxInd = 0;
		maxDeg = 0;
		for (j = 0; j < n; ++j){
			if (nodeMark[j] && nodeDegree[j] > maxDeg){
				maxDeg = nodeDegree[j];
				maxInd = j;
			}
		}
                nodeMark[maxInd] = false;
                totalCost++;
                e = hg.getNode(maxInd);
                degree[i] = degree[i-1]+nodeDegree[maxInd];
		deg += nodeDegree[maxInd];
                seeds.push_back(maxInd);
                for (j = 0; j < e.size(); ++j){
                        if (edgeMark[e[j]]){
                                continue;
                        }
                        nList = hg.getEdge(e[j]);
                        for (l = 0; l < nList.size(); ++l){
                                nodeDegree[nList[l]]--;
                        }
                        edgeMark[e[j]] = true;
                }
                i++;
        }
	for (i = 0; i < k; ++i){
                maxInd = 0;
                maxDeg = 0;
                for (j = 0; j < n; ++j){
                        if (nodeMark[j] && nodeDegree[j] > maxDeg){
                                maxDeg = nodeDegree[j];
                                maxInd = j;
                        }
                }
                nodeMark[maxInd] = false;
                ob += nodeDegree[maxInd];
	}
	ob = deg/(deg+ob);
        vector<unsigned int>().swap(nodeDegree);
        vector<int>().swap(e);
        vector<int>().swap(nList);
        vector<bool>().swap(edgeMark);
	return deg;
}

#endif
