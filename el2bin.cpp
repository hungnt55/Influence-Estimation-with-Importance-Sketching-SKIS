/*
* Functionality: convert from a graph to a binary file
* Syntax:
	./el2bin <graph file input> <binary graph output>
*/
#include <cstdio>
#include <fstream>
#include <random>
#include "GLib.hpp"
#include <cmath>
#include <cstring>

using namespace std;

int main(int argc, char ** argv)
{
	ifstream in(argv[1]);
	srand(time(NULL));
	int n,u,v;
	long long m;
	float w;
	in >> n >> m;
	printf("%d %lld\n", n, m);
	vector<int> degree(n+1,0);
	vector<vector<int> > eList(n+1);
	vector<vector<float> > weight(n+1);
	vector<float> weightR(n+1,0);

	printf("Reading the graph!\n");

	for (long long i = 0; i < m; ++i){
		in >> u >> v >> w;
		degree[v]++;
		eList[v].push_back(u);
		weight[v].push_back(w);
		weightR[u] += 1;
	}
	
	in.close();

        vector<size_t> idx(n);

	FILE * pFile;
	pFile = fopen(argv[2],"wb");
	fwrite(&n, sizeof(int), 1, pFile);
	fwrite(&m, sizeof(long long), 1, pFile);

        for (int i = 0; i < n; ++i){
		idx[i] = i;
	}
	vector<int> inv_idx(n);
	for (int i = 0; i < n; ++i){
		inv_idx[idx[i]]	= i;
	}
	
	vector<int> iTmp(n);
	
	for (int i = 0; i < n; ++i){
		iTmp[i] = degree[idx[i]+1];
	}
	
	// Write node degrees
	fwrite(&iTmp[0], sizeof(int), n, pFile);
	
	for (int i = 1; i <= n; ++i){
		// Write neighbors
		for (unsigned int j = 0; j < eList[idx[i-1]+1].size(); ++j){
			iTmp[j] = inv_idx[eList[idx[i-1]+1][j]-1]+1;
		}
		fwrite(&iTmp[0], sizeof(int), eList[idx[i-1]+1].size(), pFile);
	}

	for (int i = 1; i <= n; ++i){
		// Write weights
                fwrite(&weight[idx[i-1] + 1][0], sizeof(float), weight[idx[i-1]+1].size(), pFile);
        }

	fclose(pFile);
	printf("Done!\n");
	return 1;
}
