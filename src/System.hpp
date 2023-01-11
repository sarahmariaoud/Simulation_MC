#ifndef system_h
#define system_h

#include <iostream>
#include <array>
#include <vector>

#include "include/LowerTriangle.hpp"
#include "include/RandomObject.hpp"
#include "include/utils.hpp"


//-------------------------------------------------------------------------------
//					MULTICOMPONENT SYSTEM
//-------------------------------------------------------------------------------
class System{
private:
//---misc---//
	// RandomObject *ro;
	RandomObject *ro; 	//pointer to created random object (could possibly intilize int the system??)
//----Hyperparameters (input) ----//
	int N_types;
	std::vector<int> N;
	std::vector<std::vector<double>> F;
	std::vector<double> q;

//----system state---//
	std::vector<std::vector<int>> cs;
//--variables to keep track of----//
	double time;
	std::vector<int> N_monomers;
	int N_total;

//----interaction related stuff ----//
	std::vector<std::vector<int>> is_monomer;
	LowerTriangle<double> cp; //coalesence propensity
//-----------------------FUNCTIONS-----------------------------//
public:
//initializations
	System(int Ni,int N_types_,int *N_, double **F_, double* q_,RandomObject &ro_ );
	void initializeSystem();
//gillespie stuff
	bool do_interaction();
//OTHER GETTERS
	double get_time();
	int get_N_types();
	int get_N_total();
	int get_N_clusters();
	int get_N(int type);
	int get_N_m(int type);
	std::vector<std::vector<int>>* get_clusters();



//printing functions:
	void print_hyper();
	void print_clusters();
	void print_cp();

private:
//interactions
	void aggregate(int index_c1, int index_c2);
	void add_monomer(int type);
	void remove_monomer(int type);

//rate calculation stuff
	double coal_propensity(int c1,int c2); 
	double R_agg();
	double mon_in_propensity(int type);
	double mon_out_propensity(int type);
};



#endif