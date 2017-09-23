PARA = -std=c++11 -Wall -O3 -g -pg 
com: rwgraph.cpp skinfest.cpp el2bin.cpp option.cpp GLib.hpp
	g++ -c GLib.hpp -o GLib.o $(PARA)
	g++ -c mappedHeap.hpp -o mappedHeap.o $(PARA)
	g++ -c HeapData.hpp -o HeadData.o $(PARA)
	g++ -c option.cpp -o option.o $(PARA)
	g++ -c rwgraph.cpp -o rwgraph.o -fopenmp $(PARA)
	g++ skinfest.cpp rwgraph.o option.o -o skinfest -fopenmp $(PARA) sfmt/SFMT.c
	g++ el2bin.cpp -o el2bin $(PARA)
