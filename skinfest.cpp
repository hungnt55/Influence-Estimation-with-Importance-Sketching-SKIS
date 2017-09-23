#include "option.h"
#include "hypergraph.hpp"
#include "getMem.hpp"
#include "xorshift.hpp"
#include "sfmt/SFMT.h"
#include <iostream>
#include <sys/time.h>
#include <cmath>
#include <cstring>

using namespace std;

int main(int argc, char ** argv)
{
	srand(time(NULL));

	OptionParser op(argc, argv);
	if (!op.validCheck()){
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}
	char * inFile = op.getPara("-i");
	if (inFile == NULL){
		inFile = (char*)"network.bin";
	}

        char * model = op.getPara("-m");
        if (model == NULL)
                model = (char *) "LT";

	float scale = 1;
	char * scaledown = op.getPara("-sd");
	if (scaledown != NULL){
		scale = atof(scaledown);
	}

	Graph g;
	if (strcmp(model, "LT") == 0){
		g.readGraphLT(inFile,scale);
	} else if (strcmp(model, "IC") == 0){
		g.readGraphIC(inFile);
	} else {
		printf("Incorrect model option!");
		return -1;
	}

	int n = g.getSize();
	int m = g.getEdge();

	char * tmp;
/*	tmp = op.getPara("-epsilon");
	float epsilon = 0.1;
	if (tmp != NULL){
		epsilon = atof(tmp);
	}

	float delta = 1.0/n;
	
	tmp = op.getPara("-delta");
	if (tmp != NULL){
		delta = atof(tmp);
	}
*/
	int t = 1;
	tmp = op.getPara("-t");
	if (tmp != NULL){
		t = atoi(tmp);
	}

	HyperGraph hg(n,m);

	vector<sfmt_t> sfmtSeed = vector<sfmt_t>(t+1);
        for (int i = 0; i <= t; ++i){
                sfmt_init_gen_rand(&sfmtSeed[i], rand());
        }

	timespec start,stop;

	double h = 1;
	tmp = op.getPara("-h");
	if (tmp != NULL){
		h = atof(tmp);
	}

	long long int R = (long long int) h*n*log(n);

	long long int count = 0;

	double totalProb = g.getTotalProb();
	// cout << "Total probability: " << totalProb << endl;
	int premem = getMemValue();

        clock_gettime(CLOCK_REALTIME, &start);
	omp_set_num_threads(t);
        #pragma omp parallel
        {
                int id = omp_get_thread_num();
	        // cout << "Thread: " << id << endl;
		vector<bool> visit(n+1,false);
	        vector<int> visit_mark(n,0);
	        unsigned int num_marked = 0;

		while (count < R){
			num_marked = 0;
        		hg.pollingIC(g,visit,visit_mark,num_marked,sfmtSeed[id]);
			#pragma omp atomic
			count += num_marked;
		}
        }
	hg.updateDeg();
        clock_gettime(CLOCK_REALTIME, &stop);
        cout << "Indexing time: " << ((stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec)/exp(9*log(10))) << endl;
	cout << "Index Memory: " << (getMemValue()-premem)/1024.0 << endl;
	cout << "Total Memory: " << getMemValue()/1024.0 << endl << endl;

	// answering queries
	vector<int> seeds;
	char * list = op.getPara("-l");
	if (list == NULL){
		cout << "No list!" << endl;
		return 1;
	}
	ifstream listin(list);
	char line[256];
	unsigned int totalTime = 0;
	unsigned int estimationCount = 0;
	R = hg.getNumEdge();
	while (true){
		// read the seed set
		seeds.clear();
		listin.getline(line,256);
		while (line != NULL && strcmp(line,"") != 0 && strcmp(line," ") != 0){
			seeds.push_back(atoi(line));
			listin.getline(line,256);
		}
		if (seeds.size() <= 0){
			break;
		}
		clock_gettime(CLOCK_REALTIME, &start);
		cout << "Influence: " << (hg.estimateInf(seeds)/(double)R)*totalProb + seeds.size() << endl;
		clock_gettime(CLOCK_REALTIME, &stop);
		cout << "Time: " << ((stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec)/exp(6*log(10))) << endl << endl;
		totalTime += (stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec);
		estimationCount++;
	}
	
	cout << "Everage time: " << totalTime/exp(6*log(10))/estimationCount << endl;

	listin.close();
}
