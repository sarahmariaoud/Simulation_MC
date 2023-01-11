//	File to write data to binary
//	possible extentions:
//		-make the option for csv readable files


#ifndef file_system_h
#define file_system_h

#include <iostream>
#include <map>
#include <math.h>
#include <numeric>
#include <vector>
#include <fstream>

#include "System.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////
//									GLOBAL VARIABLES
/////////////////////////////////////////////////////////////////////////////////////////////////////

#define N_STEPS_MAX  100000

/////////////////////////////////////////////////////////////////////////////////////////////////////
//						  Binary file writting structure
//////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------
// The structure of the macro file is as follows:
//-------------------------------------------------
//	Ntypes (int_32,1)
//	('time',double,(1,)) ('N_total',int,(1,)) ('Nc',int,(1,)) ('N',int,(Ntypes,)) ('N_m',int,(Ntypes,))
//  ....
// 	....
//	('time',double,(1,)) ('N_total',int,(1,)) ('Nc',int,(1,)) ('N',int,(Ntypes,)) ('N_m',int,(Ntypes,))
//
//-------------------------------------------------
// The structure of the micro file is as follows:
//-------------------------------------------------
// ('Ntypes', int, (1,))
// ('time',double,(1,)) (C,int,(Ntypes,))
// ....
// ....
//('time',double,(1,)) (C,int,(Ntypes,))
//////////////////////////////////////////////////////////////////////////////////////////////////////
class DF{
private:
	int N_types;
	std::string fileName; //file name
	std::string fileName_micro ;
	std::string fileName_macro;

	bool write_micro;
	std::ofstream ofp_macro;
	std::ofstream ofp_micro;

public:
	DF(int N_types_,std::string &fileName_, bool write_micro_);
	~DF();
	void write(System &s);
	void writeMacro(System &s);
	void writeMicro(System &s);


};
//file name input should be without extension
inline DF::DF(int N_types_,std::string &fileName_, bool write_micro_): fileName(fileName_),write_micro(write_micro_){

 	N_types = N_types_;
	fileName_micro = fileName+".micro.dat";
	fileName_macro = fileName+".macro.dat";

	int temp_ntypes= N_types;
	ofp_macro.open(fileName_macro,std::ios::out | std::ios::binary);
	ofp_macro.write((char*)&temp_ntypes,sizeof(int));

	if(write_micro){
		ofp_micro.open(fileName_micro,std::ios::out | std::ios::binary);
		ofp_micro.write((char*)&temp_ntypes,sizeof(int));
	}

}
inline DF::~DF(){
	ofp_macro.close();
	if(write_micro){
		ofp_micro.close();
	}
}
inline void DF::write(System &s){
	writeMacro(s);
	if(write_micro) writeMicro(s);

}
inline void DF::writeMacro(System &s){
	double time = s.get_time();
	int size = 2+(2*N_types);
	int Ns[size];
	Ns[0] = s.get_N_total();
	Ns[1] = s.get_N_clusters();
	int index1 = 2, index2 = index1+N_types;
	for(int i=0; i<N_types; i++){
		Ns[index1+i] = s.get_N(i);
		Ns[index2+i] = s.get_N_m(i);
	}
	ofp_macro.write((char*)&time,sizeof(double));
	ofp_macro.write((char*)&Ns,size*sizeof(int));
}

//TODO can possibly use a diff method to write reducing number of bytes
//		ntypes ndata [time (ndata,)] [cs.size(),(ndata)]
//		[C1 (ntypes)] .. [C_cssize(0),(ntypes)]
//		[C1].....[C_cssize(1),(ntypes)]
//.....
//
void DF::writeMicro(System &s){
	
	double time = s.get_time();
	std::vector<std::vector<int>> cs = *s.get_clusters();
	for(int i=0; i<cs.size(); i++){
		std::vector<int> temp = cs[i];
		ofp_micro.write((char*)&time,sizeof(double));
		ofp_micro.write((char*)&temp[0], N_types*sizeof(int));
	}
}




#endif