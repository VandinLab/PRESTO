#include "Snap.h"

#include "temporalmotifsamplerparallel.h"
#include <random>
#include "ctpl.h"

#include <vector>
#include <algorithm>
#include <iomanip>



///////////////////////////////////////////////////////////////////////////////
TempMotifSamplerParallel::TempMotifSamplerParallel(std::string filenameG, std::string filenameM, bool parallel)
{
	if(!parallel)
	{
		std::ifstream fileG(filenameG);
		const char* DELIM = " ";
		char currentString[150]; 
		if(fileG.is_open())
		{
			long long ID = 0;
			int self_edges = 0;
			while(fileG.getline(currentString, 150, '\n'))
			{
				char* chunk = strtok(currentString, DELIM);
				int src = atoi(chunk);
				chunk = strtok(NULL, DELIM);
				int dst = atoi(chunk);
				chunk = strtok(NULL, DELIM);
				long long int timestamp = atoll(chunk);
				chunk = strtok(NULL, DELIM);
				if(src!=dst)
				{
					TEdge edge = {.src=src, .dst=dst, .tim=timestamp, .id=ID++};
					edges_.push_back(edge);
				}
				else
					self_edges++;
				delete chunk;
			}
		}
		fileG.close();
		
		std::sort(edges_.begin(), edges_.end());
		std::ifstream fileM(filenameM);
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
					std::cerr << "ERROR, motif is not in format src dst order" << std::endl;
					exit(1);
				}
				else
				{
					// Assume the edges are already ordered
					int src =stoi(results[0]);
					int dst =stoi(results[1]);
					if(src != dst)
					{
						edgesM_.push_back(std::make_pair(src, dst));
						l_++; // counting the number of edges
						max_nodes = std::max(max_nodes, src + 1);
						max_nodes = std::max(max_nodes, dst + 1);
					}
				}
			}
			Vm_ = max_nodes;
		}

		fileM.close();

	}
}

double TempMotifSamplerParallel::get_wall_time()
{
	struct timeval time;
	if (gettimeofday(&time,NULL)){
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


double TempMotifSamplerParallel::ApproximateCountMotifsSlidingWindowSkipSingleThread(int delta, double c, double mult, int b, int* samples)
{

        const double window = c * delta;
		std::cout << " c is: " << c << ", delta is: " << delta << ", r is: " << mult << '\n';
		std::random_device rd;
		std::mt19937 generator(rd());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        std::uniform_int_distribution<long long int> distributions(0, static_cast<int>(window));
        
        const int n_trials = b;

        double tot_estimate = 0;
		int times = 0;
		std::vector<std::tuple<int*,int>> data;
        for (int k = 0; k < n_trials; k++) {
            long long int offset = rand() % ((int) (window));

            double res = 0;
            long long int nxt = edges_[0].tim + offset;
		   	unsigned int idxL = 0, idxR = 0;
            while (idxR < edges_.size() && edges_[idxR++].tim < nxt);
            while (idxL < edges_.size()) {
                int nedges = idxR - idxL;

                double pi = std::min(mult * double(nedges) / edges_.size(), 1.0);
                double p = distribution(generator);
                if (p <= pi) {
                    times++;

					std::tuple<unsigned int, unsigned int, unsigned int> lims{0, 0, 0};
					
					res += ProcessSample(k,delta, window, 1, 0, idxL, idxR,std::cref(lims),
										Vm_, edgesM_, cref(edges_)) / pi;
                }

                idxL = idxR;
                if (idxR < edges_.size()) {
                    long long int t = edges_[idxR].tim;
                    nxt = t + static_cast<long long>(window);
                }
                while (idxR < edges_.size() && edges_[idxR++].tim < nxt);
            }
			tot_estimate += res;
        }
		*samples = times;
		double result = tot_estimate / double(n_trials);
        return result;
}

double TempMotifSamplerParallel::PrestoA(int delta, double c, int scale, int threads)
{
		long long delta2 = 0;
		ComputeDelta(&delta2, c, delta);
		if(delta2 == 0) {return -1.0;}
		int s = 0;
		s = scale; 

		long long min_t =  floor(edges_[l_-1].tim - (c*delta));
		long long max_t = edges_[edges_.size()-1-l_].tim; 


		// Set the random generator
		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_int_distribution<long long> randomTime(min_t, max_t);
		

		double approx = 0;
		// Computing the estimate
        std::vector<std::future<double>> futures;
        std::vector<double> probs;
		ctpl::thread_pool p(threads);

		for(int i = 0; i < s; i++)
		{
			long long t_lower = randomTime(generator);
			long long t_upper = ceil(t_lower + (c * delta));
			// TODO: use a hashmap to speedup the computation of the sample
			// Computing the i-th sample
			long long bound = ceil(t_lower);
			TEdge target = {0,0, bound, 0};
			long long b1 = ceil(t_upper);
			TEdge targ = {0, 0, b1, 0};
			unsigned int first = lower_bound(edges_.begin(), edges_.end(), target) - edges_.begin();
			unsigned int last = upper_bound(edges_.begin()+first, edges_.end(), targ) - edges_.begin();
			std::tuple<unsigned int, unsigned int, unsigned int> lims{0, 0, 0};
			futures.push_back(
							p.push(ProcessSample, delta, c, 2, delta2, first, last,std::cref(lims),
								Vm_, edgesM_, cref(edges_)));
		}

		for (int i = 0; i < (int) futures.size(); i++) {
				double part_count = futures[i].get();
                approx += part_count; 
            }
		double finalEstimate = approx/s;
		return finalEstimate;
}

double TempMotifSamplerParallel::PrestoASingThr(int delta, double c, int scale, double time)
{
		// Computing the sample size
		long long delta2 = 0;
		ComputeDelta(&delta2, c, delta);
		if(delta2 == 0) {return -1.0;}

		long long min_t =  floor(edges_[l_-1].tim - (c*delta));
		long long max_t = edges_[edges_.size()-1-l_].tim; 

		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_int_distribution<long long> randomTime(min_t, max_t);
		

		double approx = 0;
		// Computing the estimate
        std::vector<double> probs;

		int new_samps = 0;
		double st{ get_wall_time() };
		
		int id=0;
		while( (get_wall_time() - st) < time)
		{	
            //vector<TEdge>* tmp = new vector<TEdge>();
			//cout << randomTime(generator) << endl;
			long long t_lower = randomTime(generator);
			long long t_upper = ceil(t_lower + (c * delta));
			// TODO: use a hashmap to speedup the computation of the sample
			// Computing the i-th sample
			long long bound = ceil(t_lower);
			TEdge target = {0,0, bound, 0};
			long long b1 = ceil(t_upper);
			TEdge targ = {0, 0, b1, 0};
			unsigned int first = lower_bound(edges_.begin(), edges_.end(), target) - edges_.begin();
			unsigned int last = upper_bound(edges_.begin()+first, edges_.end(), targ) - edges_.begin();

			std::tuple<unsigned int, unsigned int, unsigned int> lims{0, 0, 0};
            approx += ProcessSample(id, delta, c, 2, delta2, first, last,std::cref(lims),
								Vm_, edgesM_, cref(edges_));

			new_samps++;
		
		}

		double finalEstimate = approx/new_samps;
		return finalEstimate;
}


double TempMotifSamplerParallel::PrestoE(int delta, double c, int scale, int threads)
{
		long long delta1 = 0;
		unsigned int lastpos = 0;
		ComputeDelta1(&delta1, &lastpos, c, delta);
		
		if(delta1 == 0) {return -1.0;}

		int s = scale;
		
		// Set the random generator
		
		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_int_distribution<unsigned int> randomEdge(0, lastpos);
		
		double approx = 0;
		// Computing the estimate
		std::vector<std::future<double>> futures;
		std::vector<double> probs;
		ctpl::thread_pool p(threads);
		std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> debug;
		for(int i = 0; i < s; i++)
		{
			unsigned int start = randomEdge(generator);	
			timestamp_t t_choice = edges_[start].tim;
			//timestamp_t t_lower = edges_[start].tim - static_cast<long long>(floor(c*delta));
			//unsigned int idx_lower = std::lower_bound(edges_.begin(), edges_.end(), t_lower) - edges_.begin();
			unsigned int idx_same_start = std::lower_bound(edges_.begin(), edges_.end(), t_choice) - edges_.begin();
			timestamp_t t_upper = edges_[start].tim + static_cast<long long>(ceil(c*delta));
			unsigned int end = upper_bound(edges_.begin(), edges_.end(), t_upper) - edges_.begin();
			std::tuple<unsigned int, unsigned int, unsigned int> lims{0,0,0};
			
			futures.push_back( 
							p.push(ProcessSample, delta, c, 3, delta1 , idx_same_start, end,std::cref(lims),
								Vm_, edgesM_, cref(edges_)));

		}
		
		for (int i = 0; i < (int) futures.size(); i++) {
				double part_count = futures[i].get();
				approx += part_count;
            }
		double finalEstimate = approx/s;
		return finalEstimate;
}

double TempMotifSamplerParallel::PrestoESingThr(int delta, double c, int scale, double time)
{
		// Computing the sample size
		long long delta1 = 0;
		unsigned int lastpos = 0;
		ComputeDelta1(&delta1, &lastpos, c, delta);
		
		if(delta1 == 0) {return -1.0;}

		// Set the random generator
		
		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_int_distribution<unsigned int> randomEdge(0, lastpos);
		
		double approx = 0;
		// Computing the estimate
		std::vector<double> probs;
		ctpl::thread_pool p(THREADS);

		double st{ get_wall_time() };
		
		int new_samp = 0;
		while( (get_wall_time() - st) < time)
		{	
			unsigned int start = randomEdge(generator);	
			// AFTER PAPER HANDLING MULTIPLE TIMESTAMPS..
			timestamp_t t_choice = edges_[start].tim;
			timestamp_t t_lower = edges_[start].tim - static_cast<long long>(floor(c*delta));
			
			timestamp_t t_upper = edges_[start].tim + static_cast<long long>(ceil(c*delta));
		
			unsigned int idx_lower = std::lower_bound(edges_.begin(), edges_.begin()+start, t_lower) - edges_.begin();
			unsigned int idx_same_start = std::lower_bound(edges_.begin()+idx_lower, edges_.begin()+start, t_choice) - edges_.begin();
			unsigned int end = upper_bound(edges_.begin()+start, edges_.end(), t_upper) - edges_.begin();
			std::tuple<unsigned int, unsigned int, unsigned int> lims{idx_lower, idx_same_start, end};
			approx += ProcessSample(1, delta, c, 7, delta1 , idx_same_start, end,std::cref(lims),
								Vm_, edgesM_, cref(edges_));
			new_samp++;

		}
		
		double finalEstimate = approx/new_samp;
		return finalEstimate;
}

void TempMotifSamplerParallel::ComputeDelta1(long long* delta1, unsigned int* lastpos, double c, int delta)
{	
		TEdge curr = {0, 0, static_cast<long long>(edges_.back().tim - (c*delta)), 0};
		unsigned int idx = std::lower_bound(edges_.begin(), edges_.end(), curr) - edges_.begin();
		*lastpos = idx;
		*delta1 = idx+1; // Recall the array goes from 0,...,|E|-1
}
void TempMotifSamplerParallel::ComputeDelta(long long* delta2, double c, int delta)
{
		if(edges_.size() - l_ < 0)
		{
			std::cerr << "Graph is too small" << std::endl;
		}
		*delta2 = floor(edges_[edges_.size()-l_].tim - edges_[l_-1].tim + (c * delta));
}
