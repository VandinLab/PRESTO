#ifndef snap_temporalmotifsexact_h
#define snap_temporalmotifsexact_h

#define VERBOSE 100 //Low verbose very few informative prints, high verbose lot of prints for debug and info purposes 

#include "Snap.h"

#include "temporalmotiftypes.h"

#include <vector>
#include <map>
#include <stack>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <random>

#include <cmath>

//MY LIB
#include <iomanip>
#include <functional>

// TODO: no global namespaces
using namespace std;

// TODO: REPLACE THIS WITH SNAP STRUCTURES
struct TEdge {
	long long src, dst, tim, id;
	const bool operator<(const TEdge& o) const {
		if (o.tim != tim) return tim < o.tim;
		if (o.src != src) return src < o.src;
		if (o.dst != dst) return dst < o.dst;
		if (o.id != id) return id < o.id;
		return false;
	}

	friend ostream& operator<<(ostream& in, const TEdge& o) {
		in << "(" << o.src << "," << o.dst << "," << o.tim << "," << o.id << ")";
		return in;
	}
};

class TempMotifsExact {
public:
    //TempMotifSamplerRig(const TStr& filenameG, const TStr& filenameM);
    TempMotifsExact(string filenameG, string filenameM);

    double ExactEnumerateMotifs(int delta, double window, vector<TEdge>& edges) {
    	int nV = adj_list_.size();
		// Cleanup of the data structures
    	edgeCount_.clear();
    	edgeCount_.resize(nV, 0);
    	mapGM_.clear();
    	mapGM_.resize(nV, -1);
    	mapMG_.clear();
    	mapMG_.resize(Vm_, -1);

    	double res = 0;
    	vector<int> eStack;

    	int eG = 0, eM = 0, uG = -1, vG = -1, uM = -1, vM = -1;
    	const long long int INF = numeric_limits<long long>::max();
    	long long int _t = INF;
    	while (true) {
    		int last = eStack.empty() ? -1 : eStack.back();
    		eG = FindNextMatch(eM, eG, mapMG_, mapGM_, _t, edges, last);
    		TEdge edge = {0, 0, 0, INF};
    		if (eG < (int) edges.size()) {
    			if (eM == int(edgesM_.size()) - 1) {
    				// Apply the weight function
    				if (window < 0) 
					{
						res += 1;
    				}
					else 
					{
						// Fixed Bug
						res += 1. / (1 - double(edges[eG].tim - edges[eStack[0]].tim) / window);
					}
    			} else {
    				edge = edges[eG];
    				uG = edge.src, vG = edge.dst;
    				uM = edgesM_[eM].first, vM = edgesM_[eM].second;

					mapGM_[uG] = uM;
    				mapGM_[vG] = vM;
    				mapMG_[uM] = uG;
    				mapMG_[vM] = vG;

    				edgeCount_[uG] += 1;
    				edgeCount_[vG] += 1;
    				if (eStack.empty()) {
    					_t = edge.tim + delta;
    				}
    				eStack.push_back(eG);
    				eM += 1;
    			}
    		}
    		eG += 1;

    		while (eG >= (int) edges.size() || edge.tim > _t) {
    			if (!eStack.empty()) {
    				eG = eStack.back() + 1;
    				eStack.pop_back();

    				edge = edges[eG - 1];
    				uG = edge.src, vG = edge.dst;
    				uM = edgesM_[eM].first, vM = edgesM_[eM].second;

    				if (eStack.empty()) {
    					_t = INF;
    				}
    				edgeCount_[uG] -= 1;
    				edgeCount_[vG] -= 1;

    				if (edgeCount_[uG] == 0) {
    					uM = mapGM_[uG];
    					mapMG_[uM] = -1;
    					mapGM_[uG] = -1;
    				}

    				if (edgeCount_[vG] == 0) {
    					vM = mapGM_[vG];
    					mapMG_[vM] = -1;
    					mapGM_[vG] = -1;
    				}

    				eM -= 1;
    			} else {
    				return res;
    			}
    		}
    	}
    }

    inline int FindNextMatch(
    	int eM, 
    	int eG, 
    	vector<int>& mapMG, 
    	vector<int>& mapGM, 
    	long long int _t,
    	vector<TEdge>& edges,
    	int last) {

    	int uM, vM, uG, vG;
    	uM = edgesM_[eM].first;
    	vM = edgesM_[eM].second;
    	uG = mapMG_[uM];
    	vG = mapMG_[vM];

    	vector<TEdge>* S;
    	int head = 0;
    	if (uG >= 0 && vG >= 0) {
		// Both nodes were matching: pointing to the edges in common
		// (uG, vG, t) ...
    		S = &adjMap_[uG][vG];
    	} else if (uG >= 0) {
		// Only uG was assigned to match a node in the motif
		// pointing to all the edges (uG, X, t)
    		S = &adj_list_[uG];
    	} else if (vG >= 0) {
		// Only vG was assigned to match a node in the motif
		// pointing to all the edges (X, vG, t)
    		S = &revadj_list_[vG];
    		//S = &adj_list_[vG]; 
    	} else {
    		S = &edges;
    	}
    	bool small = (S->size() < 16);
    	if (!small) {
	    	head = lower_bound(S->begin(), S->end(), edges[eG]) - S->begin();
	    }

    	for (int i = head; i < (int) S->size(); i++) {
    		auto edge = (*S)[i];
    		if (edge.id < eG || edge.tim > _t) {
    			if (small) continue;
    			else break;
    		}
			// This case should not happend if edges are sorted correctly
    		if (last != -1 && edge.tim <= edges[last].tim) {
    			continue;
    		}

    		int _eG = edge.id;
    		int _uG = edge.src, _vG = edge.dst;
    		if (uG == _uG || ((uG < 0) && (mapGM_[_uG] < 0))) {
    			if (vG == _vG || ((vG < 0) && (mapGM_[_vG] < 0))) {
    				return _eG;
    			}
    		} 
    	}
		return edges.size();
    }

    double ExactCountMotifs(int delta) {
    	InitializeStructures(edges_);
    	return ExactEnumerateMotifs(delta, -1, edges_);
    }

private:
	void InitializeStructures(vector< TEdge >& edgeList) {
		// Remapping all the nodes forming the edges
		// from 0 to n-1-----------------------
		unordered_map<int, int> remap;
		int id = 0;
		for (auto e : edgeList) {
			if (!remap.count(e.src)) {
				remap[e.src] = id++;
			}
			if (!remap.count(e.dst)) {
				remap[e.dst] = id++;
			}
		}
		for (auto& e : edgeList) {
			e.src = remap[e.src];
			e.dst = remap[e.dst];
		}
		//--------------------------------------
		// Sorting by increaing timestamps the edges
		sort(edgeList.begin(), edgeList.end());

		// Assigning an id to edges according to the result of the sorting
		for (int i = 0; i < (int) edgeList.size(); i++) {
			edgeList[i].id = i;
		}

		// TODO: IMPLEMENT SAMPLING MORE ELEGANTLY
		// Clean up of data structures
		adj_list_.clear();
		revadj_list_.clear();
		adjMap_.clear();

		// Setting the size to |E| and populating the structures
		// according to their definition
		adj_list_.resize(id);
		revadj_list_.resize(id);
		adjMap_.resize(id);
		for (const auto& edge : edgeList) {
			adj_list_[edge.src].push_back(edge);
			revadj_list_[edge.dst].push_back(edge);
			adjMap_[edge.src][edge.dst].push_back(edge);
		}
		
	} 
	
	
	//--------Edge index based adjacency list-----------
	// For each node it mantains the edges starting from "edge_i:(i,j,t)"
	vector< vector< TEdge > > adj_list_;
	// For each node it mantains the edges coming to "edge_j:(i,j,t)"
	vector< vector< TEdge > > revadj_list_;
	// vector<map_i<j, (i,j,t)>> where (i,j,t) is an edge
	vector< unordered_map<int, vector<TEdge> > > adjMap_;
	// for each node u keep the list of indexes of (u,v,t) or (v,u,t)
	vector< vector<int> > node_edge_index_;

	//-----------Structures for algo-------------------
	// maps every node of G -> # of edges assigned to them during the alg.
	vector<int> edgeCount_;
	// According to P.Mackey et al. this is an array of length |V_G|
	// of mapping nodes in G -> M
	vector<int> mapGM_;
	// According to P.Mackey et al. this is an array of length |V_M|
	// of mapping nodes in M -> G
	vector<int> mapMG_;

	// TODO: REPLACE WITH SNAP STRUCTURES
	vector< TEdge > edges_;
	// Number of nodes in the motif
	int Vm_;
	// Representation of the motif
	vector< pair<int, int> > edgesM_;
	
	//Number of edges in the motif
	int l_ = 0;
};

#endif  // snap_temporalmotifsexact_h
