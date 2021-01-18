#ifndef snap_temporalmotifsamplerparallel_h
#define snap_temporalmotifsamplerparallel_h

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
#include <thread>
#include <future>
#include <fstream>
#include <sstream>
#include <iterator>
#include <chrono>



#define VERBOSE 100
#define THREADS 8 
#define SINGLETHREAD 1 
// TODO: no global namespaces

//using namespace std;
using timestamp_t = long long int;

class Semaphore
{
	private:
		std::mutex m;
		std::condition_variable cv;
		//unsigned int count;
		std::atomic_int count;
	public:
		Semaphore(int n) : count(n) {}
		void notify()
		{
			std::unique_lock<std::mutex> l(m);
			count++;
			cv.notify_one();
			// the unique lock is destroyed as by documentation releasing the lock
		}
		void wait()
		{
			std::unique_lock<std::mutex> l(m);
			cv.wait(l, [this]{return count > 0;});
			count--;
		}
};

class Semaphore_waiter_notifier
{
		Semaphore &s;
	public:
		Semaphore_waiter_notifier(Semaphore &s) : s{s} { s.wait(); }
		~Semaphore_waiter_notifier() { s.notify(); }
};




template<typename R>
  bool is_ready(std::future<R> const& f)
  { return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready; }

// TODO: REPLACE THIS WITH SNAP STRUCTURES
struct TEdge {
	long long src;
   	long long dst;
	long long tim;
   	long long id;
	const bool operator<(const TEdge& o) const {
		if (o.tim != this->tim) return this->tim < o.tim;
		if (o.id != this->id) return this->id < o.id;
		return false; 
	}
	friend bool operator<(const TEdge &e, const timestamp_t &time) 
	{
		return (e.tim < time);
	}
	friend bool operator<(const timestamp_t &time ,const TEdge &e) 
	{
		return (time < e.tim);
	}

	friend std::ostream& operator<<(std::ostream& in, const TEdge& o) {
		in << "(" << o.src << "," << o.dst << "," << o.tim << "," << o.id << ")";
		return in;
	}
};

struct Partition
{
	unsigned int start;
	long tStart;
	unsigned int end;
	long tEnd;
	bool patch;

	const bool operator<(const Partition& p) const
	{
		if (tStart != p.tStart)
			return tStart < p.tStart;
		if ((tStart == p.tStart) && (tEnd != p.tEnd))
			return tEnd < p.tEnd;
		return false;
	}
};

static Semaphore thread_limiter(THREADS);
static Semaphore thread_limiter_exact(THREADS);


class TempMotifSamplerParallel {
public:
	double get_wall_time();
		std::mutex mtx;           // mutex for critical section
	struct Datum {
		// Vertex based data structures
		std::vector<int> edgeCount;
    	std::vector<int> mapGM;
    	std::vector<int> mapMG;
    	// Edge index based adjacency list
		std::vector< std::vector< TEdge > > adj_list;
		std::vector< std::vector< TEdge > > revadj_list;
		std::vector< std::unordered_map<int, std::vector<TEdge> > > adjMap;
		// for each node u keep the list of indexes of (u,v,t) or (v,u,t)
		std::vector< std::vector<int> > node_edge_index_;

		void InitializeStructures(std::vector< TEdge >& edgeList, int Vm) {
			std::unordered_map<int, int> remap;
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
			sort(edgeList.begin(), edgeList.end());

			for (int i = 0; i < (int) edgeList.size(); i++) {
				edgeList[i].id = i;
			}

			adj_list.resize(id);
			revadj_list.resize(id);
			adjMap.resize(id);
			for (const auto& edge : edgeList) {
				adj_list[edge.src].push_back(edge);
				revadj_list[edge.dst].push_back(edge);
				adjMap[edge.src][edge.dst].push_back(edge);
			}
			edgeCount.resize(id, 0);
			mapGM.resize(id, -1);
			mapMG.resize(Vm, -1);
			
		} 
	};

	TempMotifSamplerParallel(std::string filenameG, std::string filenameM, bool parallel=false);

	static double ExactEnumerateMotifsNew(
    	const int delta, 
    	const double window, 
    	std::vector<TEdge>* edgesP,
    	const int Vm,
    	std::vector< std::pair<int, int> > edgesM,
	   	const int type,
	   	const long long deltaT, 
		const std::vector< TEdge >& edges_,
	   	const double probability,
	   	std::vector<int>& occurr, const std::vector<Partition>& parts, const std::tuple<unsigned int, unsigned int, unsigned int> lims) {
		
    	std::vector<TEdge>& edges = *edgesP;
    	Datum data;
    	data.InitializeStructures(edges, Vm);

    	double res = 0;
    	std::vector<int> eStack;

    	int eG = 0, eM = 0, uG = -1, vG = -1, uM = -1, vM = -1;
    	const long long int INF = std::numeric_limits<long long>::max();
    	long long int _t = INF;
	   	while (true) {
    		int last = eStack.empty() ? -1 : eStack.back();
    		eG = FindNextMatch(eM, eG, data, _t, edges, last, edgesM);
    		TEdge edge = {0, 0, 0, INF};
    		if (eG < (int) edges.size()) {
    			if (eM == int(edgesM.size()) - 1) {
    				// Apply the weight function
    				if (type == 0) res += 1;
    				else if(type == 1)
					{
						res += 1. / (1 - double(edges[eG].tim - edges[eStack[0]].tim) / window);
					}
					else if(type == 2)
					{
						double r_u = 1. * ((window*delta) - (edges[eG].tim - edges[eStack[0]].tim));
						res += (deltaT/r_u);
					}
					else if(type == 3)
					{
						timestamp_t targettime = edges[eG].tim - static_cast<long long>(floor(window*delta));
						
						long long idx_lower = lower_bound(edges_.begin(), edges_.end(), targettime) - edges_.begin();

						timestamp_t time_upper = edges[eStack[0]].tim;
						long long idx_upper_tmp = upper_bound(edges_.begin(), edges_.end(), time_upper) - edges_.begin();

						long long r_u = idx_upper_tmp - idx_lower +1;
						res += (deltaT/(1.*static_cast<double>(r_u)));

					}
					else if(type == 7)
					{
						unsigned int l3 = std::get<0>(lims);; // lower index, mid index and upper index for case 3;
						unsigned int m3 = std::get<1>(lims);; // lower index, mid index and upper index for case 3;
						unsigned int u3 = std::get<2>(lims);; // lower index, mid index and upper index for case 3;
						timestamp_t targettime = edges[eG].tim - static_cast<long long>(floor(window*delta));
						long long idx_lower_tmp = lower_bound(edges_.begin()+l3, edges_.begin()+u3, targettime) - edges_.begin();
						timestamp_t time_upper = edges[eStack[0]].tim;
						long long idx_upper_tmp = upper_bound(edges_.begin()+m3, edges_.begin()+u3, time_upper) - edges_.begin();
						long long r_u = idx_upper_tmp - idx_lower_tmp;
						res += (deltaT/(static_cast<double>(r_u)));
					}
    			} else {
    				edge = edges[eG];
    				uG = edge.src, vG = edge.dst;
    				uM = edgesM[eM].first, vM = edgesM[eM].second;

    				data.mapGM[uG] = uM;
    				data.mapGM[vG] = vM;
    				data.mapMG[uM] = uG;
    				data.mapMG[vM] = vG;

    				data.edgeCount[uG] += 1;
    				data.edgeCount[vG] += 1;
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
    				uM = edgesM[eM].first, vM = edgesM[eM].second;

    				if (eStack.empty()) {
    					_t = INF;
    				}
    				data.edgeCount[uG] -= 1;
    				data.edgeCount[vG] -= 1;

    				if (data.edgeCount[uG] == 0) {
    					uM = data.mapGM[uG];
    					data.mapMG[uM] = -1;
    					data.mapGM[uG] = -1;
    				}

    				if (data.edgeCount[vG] == 0) {
    					vM = data.mapGM[vG];
    					data.mapMG[vM] = -1;
    					data.mapGM[vG] = -1;
    				}

    				eM -= 1;
    			} else {
    				delete edgesP;
					/*
					try{
						delete &edges;
					} catch(...){cerr << "c'è stato un errore" << endl;}
					*/
					//cerr << "delted pointer, returning partial count.. " << endl;
    				return res;
    			}
    		}
    	}
    }

    static double ExactEnumerateMotifs(
    	const int delta, 
    	const double window, 
    	std::vector<TEdge>* edgesP,
    	const int Vm,
    	std::vector< std::pair<int, int> > edgesM,
	   	const int type,
	   	const long long deltaT, 
		const std::vector< TEdge >& edges_,
	   	const double probability,
	   	std::vector<int>& occurr,
		const std::tuple<unsigned int, unsigned int, unsigned int> lims) {
	std::vector<Partition> fill;
	return ExactEnumerateMotifsNew(delta, window, edgesP, Vm, edgesM, type, deltaT, std::cref(edges_), probability, std::ref(occurr), std::cref(fill), lims);
	}


    static inline int FindNextMatch(
    	int eM, 
    	int eG, 
    	Datum& data, 
    	long long int _t,
    	std::vector<TEdge>& edges,
    	int last,
    	std::vector< std::pair<int, int> >& edgesM) {

    	int uM, vM, uG, vG;
    	uM = edgesM[eM].first;
    	vM = edgesM[eM].second;
    	uG = data.mapMG[uM];
    	vG = data.mapMG[vM];

    	std::vector<TEdge>* S;
    	int head = 0;
    	if (uG >= 0 && vG >= 0) {
    		S = &data.adjMap[uG][vG];
    	} else if (uG >= 0) {
    		S = &data.adj_list[uG];
    	} else if (vG >= 0) {
    		S = &data.revadj_list[vG];
    	} else {
    		S = &edges;
    	}
    	bool small = (S->size() < 16);
    	if (!small) {
	    	head = std::lower_bound(S->begin(), S->end(), edges[eG]) - S->begin();
	    }

    	for (int i = head; i < (int) S->size(); i++) {
    		auto edge = (*S)[i];
    		if (edge.id < eG || edge.tim > _t) {
    			if (small) continue;
    			else break;
    		}
    		if (last != -1 && edge.tim <= edges[last].tim) {
    			continue;
    		}

    		int _eG = edge.id;
    		int _uG = edge.src, _vG = edge.dst;
    		if (uG == _uG || (uG < 0 && data.mapGM[_uG] < 0)) {
    			if (vG == _vG || (vG < 0 && data.mapGM[_vG] < 0)) {
    				return _eG;
    			}
    		} 
    	}
    	return edges.size();
    }

    double ApproximateCountMotifsSlidingWindowSkipSingleThread(int delta, double c = 3e1, double mult = 3e1, int b=5, int* samples = nullptr);
   	
	static double ProcessSample(int id, double delta, double c, int type, double deltaT,unsigned int first,unsigned int last, const std::tuple<unsigned int, unsigned int, unsigned int>& lims ,
							int Vm_, std::vector<std::pair<int,int>> edgesM, const std::vector< TEdge >& edges_)
	{
        std::vector<TEdge>* tmp = new std::vector<TEdge>();
		tmp->insert(tmp->begin(), edges_.begin()+first, edges_.begin()+last); // Sample goes from [first,last)
		std::vector<int> a;	
		double result = ExactEnumerateMotifs(delta, c, tmp, Vm_, edgesM, type, deltaT, std::cref(edges_), 0, std::ref(a), lims);
		return result;
	}


	/* Sampling Version 2.. see thesis */
	double PrestoA(int delta, double c, int scale, int threads=8);

	double PrestoASingThr(int delta, double c, int scale, double time);
	
	double PrestoE(int delta, double c, int scale, int threads=8);

	double PrestoESingThr(int delta, double c, int scale, double time);

	void ComputeDelta1(long long* delta1, unsigned int* lastpos, double c, int delta);
	
	void ComputeDelta(long long* delta2, double c, int delta);
	
private:
	double get_wall_time_inner(){

		struct timeval time;
		if (gettimeofday(&time,NULL)){
			//  Handle error
			return 0;
		}
		return (double)time.tv_sec + (double)time.tv_usec * .000001;
	}
	std::vector< std::vector<TEdge> > all_components;
	std::vector< TEdge > edges_;
	std::vector< std::pair<int, int> > edgesM_;
    std::default_random_engine generatorLiu;
	int Vm_;
	int l_ = 0;
};

#endif  // snap_temporalmotifsamplerparallel_h
