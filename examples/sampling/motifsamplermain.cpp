#include "temporalmotifsamplerparallel.h"

#include <omp.h>

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <mutex>

double get_wall_time(){

	struct timeval time;
	if (gettimeofday(&time,NULL)){
		//	Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int main(int argc, char* argv[]) {
	const rlim_t kStackSize = RLIM_INFINITY;   // min stack size = 16 MB
	struct rlimit rl;
	int status;

	status = getrlimit(RLIMIT_STACK, &rl);
	if (status == 0)
	{
		if (rl.rlim_cur < kStackSize)
		{
			rl.rlim_cur = kStackSize;
			status = setrlimit(RLIMIT_STACK, &rl);
			if (status != 0)
			{
				fprintf(stderr, "setrlimit returned result = %d\n", status);
			}
		}
	}
	

  Env = TEnv(argc, argv, TNotify::StdNotify);

  const TStr temporal_graph_filename =
	Env.GetIfArgPrefixStr("-i:", "simple-example.txt",
			  "Input directed temporal graph file");
  const TStr motif_graph_filename =
	Env.GetIfArgPrefixStr("-m:", "simple-motif.txt",
		"Input directed motif graph file");
  const TFlt delta =
	Env.GetIfArgPrefixFlt("-delta:", 4096, "Time window delta");
  const int num_threads =
	Env.GetIfArgPrefixInt("-nt:", 4, "Number of threads for parallelization");
  const TFlt c =
	Env.GetIfArgPrefixFlt("-c:", 30, "Window size multiplier");
  const int type =
	Env.GetIfArgPrefixFlt("-type:", 1, "Which alg?");
  const int samp=
	Env.GetIfArgPrefixFlt("-samp:", 100, "Number of samples (s) used for PRESTO");
  const TFlt mult =
	Env.GetIfArgPrefixFlt("-r:", 30, "Sampling probability multiplier of LS");
  const int bLS =
	Env.GetIfArgPrefixFlt("-b:", 5, "Number of shifts of LS");
  
  double st = get_wall_time(); 
  TempMotifSamplerParallel tmc(temporal_graph_filename.CStr(), motif_graph_filename.CStr());

  Env.PrepArgs(TStr::Fmt("Temporalmotifs. build: %s, %s. Time: %s",
	   __TIME__, __DATE__, TExeTm::GetCurTm()));  
	printf("Time to load the dataset: %lfs\n\n", get_wall_time() - st);
		srand(time(nullptr));
  
  double init;
  double end;
	int samples = 0;
	
	
	// RUNNING PARALLEL VERSION IF TYPE == 0
	if(type == 0)
	{

		for(int i=0; i < 5; i++)
		{	
			if(samp == 0)
			{
				std::cout << "Samples have to be >1!!\n";
				exit(0);
			}
			samples = samp;

			init = get_wall_time();
			printf("--PRESTO-A-- approximation estimate is : %lf\n", tmc.PrestoA(delta, c, samples, num_threads));
			printf("--PRESTO-A-- run time: %lfs\n\n", get_wall_time() - init);

			init = get_wall_time();
			double a = tmc.PrestoE(delta, c, samples, num_threads);
			printf("--PRESTO-E-- run time : %lfs\n", get_wall_time() - init);
			printf("--PRESTO-E-- approximation estimate is : %lf\n\n", a); 
		}	
	}
	// RUNNING FIRST LS AND THEN OUR ALGORITHMS, see main text
	else if(type == 1)
	{

		for(int i=0; i < 5; i++)
		{	
			samples =0;
			init = get_wall_time();
			double st = init;
			printf("--LS-- Solution: %lf\n", tmc.ApproximateCountMotifsSlidingWindowSkipSingleThread(delta, c, mult, bLS, &samples));
			end = get_wall_time();
			printf("--LS-- run time: %lfs\n\n", end - init);

			init = get_wall_time();
			printf("--PRESTO-A-- approximation estimate is : %lf\n", tmc.PrestoASingThr(delta, c, samples, (end-st)));
			printf("--PRESTO-A-- run time: %lfs\n\n", get_wall_time() - init);

			init = get_wall_time();
			double a = tmc.PrestoESingThr(delta, c, samples, (end-st));
			printf("--PRESTO-E-- run time : %lfs\n", get_wall_time() - init);
			printf("--PRESTO-E-- approximation estimate is : %lf\n\n", a); 

		}	
	}
  return 0;
}
