#include "Snap.h"

#include "temporalmotifsexact.h"

#include <vector>
#include <algorithm>

#include <fstream>
#include <sstream>
#include <iterator>
///////////////////////////////////////////////////////////////////////////////

TempMotifsExact::TempMotifsExact(string filenameG, string filenameM)
{
	ifstream fileG(filenameG);
	if(fileG.is_open())
	{
		std::string line;
		long int ID = 0;
		int self_edges = 0;
		while(getline(fileG, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> results(std::istream_iterator<
						std::string>{iss},std::istream_iterator<std::string>());
			if(results.size() != 3)
			{
				std::cerr << "ERROR, dataset is not in format src dst timestamp" << endl;
				exit(1);
			}
			else
			{
				int src =stoi(results[0]);
				int dst =stoi(results[1]);
				long long int timestamp = stoll(results[2]);
				if(src != dst)
				{
					TEdge edge = {.src=src, .dst=dst, .tim=timestamp, .id=ID++};
					edges_.push_back(edge);
				}
				else
					self_edges++;
			}
		}
		cout << "removed self-edges is: " << self_edges << endl;
	}
	fileG.close();
	
  	//InitializeStructures(edges_);

	ifstream fileM(filenameM);
	if(fileM.is_open())
	{
		std::string line;
		int max_nodes = 0;
		while(getline(fileM, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> results(std::istream_iterator<
						std::string>{iss},std::istream_iterator<std::string>());
			if(results.size() != 3)
			{
				std::cerr << "ERROR, motif is not in format src dst order" << endl;
				exit(1);
			}
			else
			{
				// Assume the edges are already ordered
				int src =stoi(results[0]);
				int dst =stoi(results[1]);
				if(src != dst)
				{
      				edgesM_.push_back(make_pair(src, dst));
	  				l_++; // counting the number of edges
					max_nodes = max(max_nodes, src + 1);
					max_nodes = max(max_nodes, dst + 1);
				}
			}
		}
		Vm_ = max_nodes;
	}

	fileM.close();
}
