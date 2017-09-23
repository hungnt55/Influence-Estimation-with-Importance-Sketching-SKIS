# Influence-Estimation-with-Importance-Sketching-SKIS
Information:
--------------------------------------------------------
Version 1.0: Implementation of influence estimation algorithm based on the importance sketching SKIS under the Independent Cascade(IC) model. For more details about SKIS and the algorithms, please refer to our paper:

[H. T. Nguyen, T. P. Nguyen, NhatHai Phan, T. N. Dinh, Importance Sketching of Influence Dynamics in Billion-scale Networks, IEEE International Conference on Data Mining series (ICDM), 2017](https://arxiv.org/abs/1709.03565)



Terms of use:
--------------------------------------------------------
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.


Requirements:
--------------------------------------------------------
In order to compile all the tools, it requires GCC 4.7.2 and later.


Compile:
--------------------------------------------------------
Type `make' command to compile everything


How to use:
--------------------------------------------------------
This package offers a set of functions to use in order to estimate the influences of multiple seed sets contained in a file. A typical sequence of actions is as follows:

1. Conversion from a text format to binary file for fast graph reading
	./el2bin <input file> <output file>

    <input file>: the path to text file in edge list format: the first line contains the number of nodes n and number of edges m, each of the next m lines describes an edge following the format: <src> <dest> <weight>. Node index starts from 1.
    <output file>: the path to binary output file

2. Run the algorithm to estimate the ifnluences of seed sets
	./skinfest [Options]

    Options:

        -i <binary graph file>
            specify the path to the binary graph file (default: network.bin)

        -l <seed file path>
            path to the file containing the seed sets. Each seed node is listed in a line and between two seed sets, there is a blank line to separate them

	-h <h parameter>
            determine the number of samples, i.e. h*n*log(n). Refer to our paper for more details

	-m <model>
	    diffusion model (currently only IC is supported)



     Output format:
	The outputs are printed on standard output stream in the following order
		Indexing time: <time>
		Index Memory: <index mem>
		Total Memory: <total mem>

		Influence: <influence of first seed set>
		Time: <querying time of the first seed set>

		...

		Everage time: <average querying time over all seed sets>

Example:
--------------------------------------------------------
Estimate influences of two seed sets in DBLP networks. The seed set file (dblp.seeds) contains:
		1
		2
		3
		4

		5
		6
		7
		8

	1. Convert to binary file:
		./el2bin dblp_format.txt dblp.bin
	
	2. Run our algorithm with h=5:
		./skinfest -i dblp.bin -l dblp.seeds -h 5 -m IC

	The output:

		Indexing time: 9.14089
		Index Memory: 699.676
		Total Memory: 806.676

		Influence: 1491.42
		Time: 0.14088

		Influence: 1022.63
		Time: 0.067684

		Everage time: 0.104282
