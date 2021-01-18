#include "temporalmotifs.h"
#include"temporalmotifsexact.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <time.h>
#include <sys/time.h>

double get_wall_time();

double get_wall_time(){
	struct	timeval	time;
	if(gettimeofday(&time,NULL)){
	//	Handle	error
	return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int	main(int argc, char* argv[]) {
	//	SNAP data structure for	parsing	the	input lines
	Env	= TEnv(argc, argv, TNotify::StdNotify);

	const TStr temporal_graph_filename = Env.GetIfArgPrefixStr("-i:", "../test-graphs/email-Eu-core-temporal.txt", "Input directed temporal graph file");
	const TStr motif_graph_filename	= Env.GetIfArgPrefixStr("-m:", "../motifs/M-1-1.txt", "Input directed motif graph file");
	const TFlt delta = Env.GetIfArgPrefixFlt("-delta:", 86400, "Time window delta");

	TempMotifsExact	tmc(temporal_graph_filename.CStr(), motif_graph_filename.CStr());

	Env.PrepArgs(TStr::Fmt("Temporalmotifs.	build:	%s,	%s.	Time: %s",
	__TIME__, __DATE__, TExeTm::GetCurTm()));	

	//Measure time and run the exact routine
	double start = get_wall_time();
	double result = tmc.ExactCountMotifs(delta); 
	cout << "Running time of the exact routine is: " << get_wall_time()-start << endl;
	cout << "Number of delta instances is: " << result << endl;
	return 0;
}
