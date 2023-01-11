//clang++  -Xpreprocessor -fopenmp -std=c++17 main.cpp System.cpp -o main.out -lomp -Wall


#include <iostream>

#include <filesystem>

#if defined(_OPENMP)
   #include <omp.h>
#endif


#include "System.hpp"
#include "FileSystem.hpp"
#include "include/TickTock.hpp"
#include "include/utils.hpp"
#include "include/nlohmann/json.hpp"
using json = nlohmann::json;


void test();

void run_simulation(int argc,char **argv);


//////////////////////////////////////////////////////////////////////////////////////////////////////
//							CODE STARTS HERE
//////////////////////////////////////////////////////////////////////////////////////////////////////
// IF NO ARGUMENTS PROVIDED RUNS TEST
// Input should be of the form:
//	./main.out N_rels tmax Ntypes Ni[Ntypes] q[Ntypes] F[Ntypes][Ntypes] "dir/"
//		-N_rels: 				number of realizations [ integer > 0]
//		-tmax:					maximum time of simulation [double > 0]
//		-Ntypes:	
//		-Ni[Ntypes]:
//		-q[Ntypes]:		
//		-F[Ntypes][Ntypes]:
//		-"dir/":				directory to put data
//								if empty then saves in current working directory 
//////////////////////////////////////////////////////////////////////////////////////////////////////

#define N_STEPS_MAX  100000
int main(int argc, char **argv){


	if(argc==1){
    	std::cout << "inside test module" << std::endl;
    	// test_gillespie();
    	test();
	}
	else{
		std::cout << "--------------------------------------------------------------" << std::endl;
		std::cout << "Number of arguments =  " << argc << std::endl;
		for(int i=0; i<argc; i++){
			std::cout << argv[i] << " ";
		}
		std::cout << std::endl;
		std::cout << "--------------------------------------------------------------" << std::endl;

		//calculate the number of arguments needed:
		int D = std::stoi(argv[3]);
		int Nargs = 5+ (2*D) + (D*D);
		if(argc!=Nargs) throw std::invalid_argument("wrong number of arguments");

		std::cout << "--------------------------------------------------------------" << std::endl;
		std::cout << "\t\t START SIMULATION" << std::endl;
		std::cout << "--------------------------------------------------------------" << std::endl;
		Timer tc;
		tc.tick();

		run_simulation(argc,argv);

		tc.tock();
		std::cout << "--------------------------------------------------------------" << std::endl;
		std::cout << "\t END SIMULATION ( total run time =  " << tc.duration() << " ms  )" << std::endl; 
		std::cout << "--------------------------------------------------------------" << std::endl;

	}
	return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
//							RUNNING THE SIMULATION
//////////////////////////////////////////////////////////////////////////////////////////////////////
void run_simulation(int argc,char **argv){

	std::ifstream f(argv[argc-1]);
	json j = json::parse(f);
	f.close();

	const int N_types = std::stoi(argv[3]);

//1-CONVERTING THE ARGUMENTS
	std::string dir = j["data_file_location" ];

	int N_rels = std::stoi(argv[1]);
	double t_max = std::stod(argv[2]);
	int *N; N = (int*) malloc(N_types*sizeof(int));
	double *q; q= (double*) malloc(N_types*sizeof(double));
	double **F; F = (double**) malloc(N_types*sizeof(double*));
	for(int i=0; i<N_types; i++) F[i]= (double*) malloc(N_types*sizeof(double));

	int indexN = 4, indexQ= indexN+N_types, indexF = indexQ+N_types;
	for(int i=0; i<N_types; i++){
		N[i] = std::stoi(argv[indexN]);
		q[i] = std::stod(argv[indexQ]); 

		for(int j=0; j<N_types; j++){
			F[i][j] = std::stod(argv[indexF]);
			indexF++;
		}
		indexN++;
		indexQ++;


	}
	int Ni=0;
	for(int i=0; i<N_types; i++) Ni+=N[i];

//2-CHECKING IF THE DIRECTORY IS PROPER
	if (!((dir[dir.length()-1] == '/' )) && (dir.length()>0)) throw std::invalid_argument("directory is not valid ");
	
	std::string strN = std::string("")+tostr(N[0]);
	std::string strQ = std::string("")+tostr(q[0]);
	std::string strF = std::string("")+tostr(F[0][0]);
	for(int i=0; i<N_types; i++){

		if(i!=0){
			 strN = strN + "-" +tostr(N[i]);
			 strQ = strQ + "-" +tostr(q[i]);
		}
		for(int j=0; j<N_types; j++){
			if( !(i==0 && j==0)){
				strF = strF +"-"+tostr(F[i][j]);
			}
		}
	}
	std::string data_folder = dir+"Ni_"+strN+"_q_"+strQ+"_F_"+strF;
	std::cout << "WRITTING DATA TO : " << data_folder<< "/" <<std::endl;
	std::filesystem::create_directories(data_folder);
	std::string time_str = get_time_string();
	std::cout << "TIME START : " << time_str <<std::endl;

//3-CREATING TIME TRACKERS
	int time_rel[N_rels];
	int step_rel[N_rels];


//4-RUNNING THE REALIZATIONS (IN PARALLEL)
#if defined(_OPENMP)
   std::cout << "USING  [" << omp_get_max_threads() << "] THREADS" << std::endl;
   #pragma omp parallel for
#endif	 
    for(int rel =0; rel < N_rels; rel++){
	
	//A-SETTING THE TIMER
    	Timer tc;
    	tc.tick();
	
	//B-CREATING THE DATA FOLDER
    	std::string data_file = data_folder+"/"+time_str+"-"+std::to_string(rel);
    	DF df(N_types,data_file,true);
	
	//C- INITIALIZING THE SYSTEM
    	RandomObject r = RandomObject();
		System s(Ni,N_types,N,F,q,r);
		df.write(s);

		if(rel==0){
			std::cout << "--------------------------------------------------------------" << std::endl;
			s.print_hyper();
			std::cout << "--------------------------------------------------------------" << std::endl;

		}
	//D- RUNNING THE SIMULATION
		bool cont = true;
		for(int step=0; (step < N_STEPS_MAX && (s.get_time()<t_max) && cont ); step++){
			cont = s.do_interaction();
			df.write(s);
			step_rel[rel] =step;
			if(step == N_STEPS_MAX-1 ) std::cout << "MAXIMUM STEPS ACHIEVED ENDING SIMULATION" << std::endl;
		}
	//E-END TIMER AND PLACE IN TRACKER
		tc.tock();
		time_rel[rel]=tc.duration();
    }

//5-PRINTING OUT THE TIMES FOR EACH REALIZATION
    double avg_time = 0;
	std::cout << "times [  " << time_rel[0] << "  ";
	for(int i=1; i<N_rels; i++ ){
		avg_time+=time_rel[i];
		std::cout << time_rel[i] << "  ";
	}
	std::cout<< "] [ avg = " << avg_time/(1.0*N_rels) << " ] "<<std::endl;
	
	double avg_step= 0;
	std::cout << "steps [  " << step_rel[0] << "  ";
	for(int i=1; i<N_rels; i++ ){
		avg_step+= step_rel[i];
		std::cout << step_rel[i] << "  ";
	}
	std::cout<< "] [ avg = " << avg_step/(1.0*N_rels) << " ] "<<std::endl;


	//deallocation
	free(N);
	free(q);
	for(int i=0; i<N_types; i++ ) free(F[i]);
	free(F);

}


void test(){



}