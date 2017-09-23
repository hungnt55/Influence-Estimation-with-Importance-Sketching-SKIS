#ifndef _RWGRAPH_H
#define _RWGRAPH_H
#include <vector>
#include <random>
#include "sfmt/SFMT.h"

typedef uint32_t UI;
typedef uint64_t ULL;


// Graph class defines fundamental operators on graph
class Graph
{
	friend class HyperGraph;
	private:
		UI UI_MAX = 4294967295U;
		ULL ULL_MAX = 18446744073709551615ULL;

		// number of nodes
		unsigned int numNodes;
		// number of edges
		unsigned int numEdges;
		// adjacency list
		std::vector<std::vector<int> > adjList;
//		std::vector<std::vector<int> > adjIndex;
		std::vector<int> node_deg;
		std::vector<std::vector<UI> > weights;
		std::vector<UI> nodeProb;
		double totalProb;
	
	public:
		// get a vector of neighbours of node u
		const std::vector<int> & operator [] (int u) const;
		// return weights of neighbours of node u
		const std::vector<UI> & getWeight(int u) const;

		// get a vector of neighbours of node u
		const std::vector<int> & operator [] (int u);
		// return weights of neighbours of node u
		const std::vector<UI> & getWeight(int u);

		// get degree of node u
		int getDegree(int u) const;
		// get size of the graph
		int getSize() const;
		// get number of edges
		int getEdge() const;

		void findEdge(int ind, int & end1, int & end2);

		// read graph from a file
		void readGraphLT(const char * filename,float scale);
		// read graph from a file
		void readGraphIC(const char * filename);
		// write the graph to file
		void writeToFile(const char * filename);

		double getTotalProb();
};

class HyperGraph
{
	private:
		UI UI_MAX = 4294967295U;
                ULL ULL_MAX = 18446744073709551615ULL;

		// store the edges that a node is incident to
		std::vector<std::vector<int> > node_edge;
		// store hyperedges
		std::vector<std::vector<int> > edge_node;
		unsigned int curEdge;
		unsigned int maxDegree;
		unsigned int numNodes;
		unsigned int numEdges;
		sfmt_t sfmtSeed;
		inline int randIndex_bin(const std::vector<UI> &w, unsigned int si,sfmt_t & sfmtSeed, int ben, int end, bool normalized);
		inline int randIndex_lin(const std::vector<UI> &w, unsigned int si,sfmt_t & sfmtSeed);

		std::vector<bool> seeds;
		std::vector<int> activeNodes;
		double totalProb;
		std::vector<bool> edgeMask;
		std::vector<int> visitEdge;
		int num_masked;
		std::vector<int> node_root;

	public:
		HyperGraph(unsigned int n, unsigned int m);
		void updateDeg();
		void updateEdge();
		void addEdge(std::vector<int> & edge);
		void addEdgeD(std::vector<int> & edge);
		int getMaxDegree();
		const std::vector<int> & getEdge(int e) const;
		const std::vector<int> & getEdge(int e);
		const std::vector<int> & getNode(int n) const;
		const std::vector<int> & getNode(int n);
		// get number of edges
                int getNumEdge() const;
		void clearEdges();
		void pushHyperedges(std::vector<std::vector<int> > & hyperdges);
		
		void setSeeds(std::vector<int> &s);
		
		void pollingLT1(Graph &g, std::vector<bool> &visit, std::vector<int> &mark_visit, std::vector<bool> & link, std::vector<double> & prob);
		int pollingLT2(Graph &g, std::vector<unsigned int> & link, std::vector<bool> & interdict, bool & success, unsigned int k, std::vector<bool> &visit, std::vector<int> &mark_visit, std::vector<double> & prob, unsigned int & num_marked, sfmt_t & sfmtSeed);
		int pollingLT3(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &mark_visit, std::vector<double> & prob, unsigned int & num_marked, sfmt_t & sfmtSeed);
		bool pollingLT(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &mark_visit);
		void pollingIC1(Graph &g, std::vector<bool> &visit, std::vector<int> &visit_mark);
		int pollingIC2(Graph &g, std::vector<unsigned int> & link, unsigned int k, std::vector<bool> &visit, std::vector<int> &visit_mark, unsigned int & num_marked, sfmt_t & sfmtSeed1);
		void pollingIC(Graph &g, std::vector<bool> &visit, std::vector<int> &visit_mark, unsigned int & num_marked, sfmt_t & sfmtSeed1);
		int estimateInf(std::vector<int> & seeds);
};

float getCurrentMemoryUsage();

#endif
