
#include "System.hpp"
//////////////////////////////////////////////////////////////////////////////////
//								FUNCTIONS
///////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
//					INTITIALIZATION
//-----------------------------------------------------------------------------
System::System(int Ni,int N_types_,int *N_, double **F_, double* q_,RandomObject &ro_  ):cp(Ni){
	N_types = N_types_;
	ro = &ro_;
	//setting the correct size for the vectors 
	N.resize(N_types,0);
	F.resize(N_types);
	q.resize(N_types,0);
	N_monomers.resize(N_types,0);
	for(int i=0; i<N_types; i++) F[i].resize(N_types);
	//initializing parameters
	time=0;
	N_total = 0;
	for(int i=0; i<N_types; i++) {
		N_total += N_[i];
		N[i]= N_[i];
		N_monomers[i] = N_[i];
		q[i] = q_[i];
		for(int j=0; j<N_types; j++) F[i][j] = F_[i][j];
	}
	//we now need to initialize cp is_monomer	
	initializeSystem();
}

void System::initializeSystem(){
	//setting the size for is_monomer
	is_monomer.resize(N_types);
	for(int i=0; i<N_types; i++) is_monomer[i].resize(N_total,0);	
	//initializing the system
	int sum=N[0];
	int current_type=0;
	std::vector<int> temp(N_types,0);
	for(int i=0; i<N_total; i++){

		//adding a new cluster
		for(int j=0; j<N_types; j++){
			if(j==current_type) temp[j]=1;
			else temp[j]=0;
		}
		cs.push_back(temp);
		//setting the aggregation array
		for(int j=0; j<i; j++){
			cp.set(i,j,coal_propensity(i,j));
		}
		//setting is monomer tracker
		is_monomer[current_type][i]=1;
		//seeing if we have to move on to the next type
		if(i == sum-1){
			current_type++;
			sum += N[current_type];
		}

	}

}
//------------------------------------------------------------------------------
//				INTERACTIONS:AGGREGATION
//------------------------------------------------------------------------------
void System::aggregate(int index_c1, int index_c2){
	
	if(index_c1==index_c2) throw std::invalid_argument("cant aggregate the same cluster");

	int i1 = std::min(index_c1,index_c2);
	int i2 = std::max(index_c1,index_c2);

	//adding c2 to c1 and adjusting monomer trackers
	for(int i=0; i<N_types; i++) {
		cs[i1][i] += cs[i2][i];
		if(is_monomer[i][i1]==1){ //monomer of type i
			is_monomer[i][i1] = 0;
			N_monomers[i] --;
		}
		if(is_monomer[i][i2]==1){ //monomer of type i
			is_monomer[i][i2] = 0;
			N_monomers[i] --;
		}
	}

	//removing the old cluster
	cs.erase(cs.begin()+i2);
	cp.remove(i2);
	for(int i=0; i<N_types; i++) is_monomer[i].erase(is_monomer[i].begin()+i2);

	//updating aggregation matrix
	for(int i=0; i<cs.size(); i++){
		if(i1!=i) cp.set(i1,i,coal_propensity(i1,i));
	}

}
//------------------------------------------------------------------------------
//				INTERACTIONS: ADD MONOMER
//------------------------------------------------------------------------------
void System::add_monomer(int type){
	if(type>=N_types) throw std::invalid_argument("type doesnt exist");

	N[type]++;
	N_monomers[type]++;
	N_total = get_N_total();

	std::vector<int> new_monomer(N_types,0);
	for(int i=0; i<N_types; i++){
		if(i == type){ 
			new_monomer[i]=1;
			is_monomer[i].push_back(1);
		}
		else{
			 new_monomer[i]=0;
			 is_monomer[i].push_back(0);
		}
	}
	cs.push_back(new_monomer);
	cp.add();

	int index_monomer = cs.size()-1;
	for(int i=0; i<cs.size(); i++){
		if(i!= index_monomer) cp.set(index_monomer,i,coal_propensity(i,index_monomer));
	}
}
//------------------------------------------------------------------------------
//				INTERACTIONS : REMOVE MONOMER 
//------------------------------------------------------------------------------
void System::remove_monomer(int type){
	
	if(N_monomers[type]<=0) throw std::invalid_argument(" no monomers to remove");
	//randomly select a monomer of type 
	int ith_monomer_to_remove = ro->get_int(1,N_monomers[type]);
	int index_to_remove = search_exceeds_cum<int>(is_monomer[type],ith_monomer_to_remove);	
	//adjust counting vectors
	N_monomers[type]--;
	N[type]--;
	//remove it from all indexed vectors
	cs.erase(cs.begin()+index_to_remove);
	cp.remove(index_to_remove);
	for(int i=0; i<N_types; i++) is_monomer[i].erase(is_monomer[i].begin()+index_to_remove);	
}
//------------------------------------------------------------------------------
// 					PRINTING STUFF
//------------------------------------------------------------------------------
void System::print_hyper(){
	std::cout << "N = \t[\t" << N[0] << "\t";
	for(int i=1; i<N_types; i++){
		std::cout << N[i] <<"\t";
	}
	std::cout<<"]\n";
	std::cout << "q = \t[\t" << q[0] << "\t";
	for(int i=1; i<N_types; i++){
		std::cout << q[i] <<"\t";
	}
	std::cout<<"]\n";
	std::cout << "F = ";
	for(int i=0 ; i<N_types; i++){
		std::cout << "\n\t" << F[i][0] << "\t";
		for(int j=1; j<N_types; j++){
			std::cout<< F[i][j] << "\t";
		}
	}
	std::cout << std::endl;
}
void System::print_clusters(){	
	std::cout << " Number of objects : ";
	for(int i=0; i<N_types ; i++) std::cout << N[i] << "\t";
	std::cout << "\n Number of monomers : ";
	for(int i=0; i<N_types ; i++) std::cout << N_monomers[i] << "\t";
	std::cout << "\n Number of clusters : " << cs.size() << std::endl;
	
	for(int i=0; i<cs.size(); i++){
		std::cout << i << " --->  ( " << cs[i][0];
		for(int j=1; j<N_types; j++){
			std::cout<< " , "  << cs[i][j];
		}
		std::cout << " )\t";
		for(int j=0; j<N_types; j++){
			std::cout << is_monomer[j][i] << "\t";
		}
		std::cout<<std::endl;
	}
}
void System::print_cp(){
	std::cout << "aggregation propensities \n";
	cp.print_lower_triangle();	
}

//------------------------------------------------------------------------------
//				RATE STUFF
//------------------------------------------------------------------------------
double System::R_agg(){
	//its time dependent on the current size of the system:
	if (N_total <= 1) return 0;
	return  2.0/(1.0*N_total*(N_total-1)); //TODO confirm its correct
}
double System::coal_propensity(int c1, int c2){
	double prop1 =0;
	double prop2=0;
	for(int i=0; i<N_types; i++){
		for(int j=0; j<N_types; j++){
			prop1+= (cs[c1][i]*cs[c2][j]*F[i][j]) ;
			prop2+= (cs[c2][i]*cs[c1][j]*F[i][j]);
		}
	}
	return 0.5*(prop1+prop2); 
}
double System::mon_in_propensity(int type){
	if(q[type]>0) return q[type];
	else return 0;
}
double System::mon_out_propensity(int type ){
	if(q[type]<0 && N_monomers[type]>0) return -1*q[type];
	else return 0;
}
//-------------------------------------------------------------------------
//				GILLESPIE
//--------------------------------------------------------------------------
bool System::do_interaction(){

//0-TOTAL PROPENSITY
	//aggregation and cp
	double R = R_agg();
	double alpha_agg = R*cp.get_cum();
	//input/output propensity and alpha 
	//at this point its constant so no need to recalculate it but to keep things general
	std::vector<double> inp(N_types,0),outp(N_types,0);
	double alpha_in=0, alpha_out =0;
	for(int i=0; i<N_types ; i++){
		inp[i] = mon_in_propensity(i);
		alpha_in +=inp[i];
		outp[i] = mon_out_propensity(i);
		alpha_out += outp[i];
	}
	double alpha_tot = alpha_agg+alpha_in+alpha_out;
	if (isZero(alpha_tot)){
		return false; //if the propensities reach zero then we want to end the simulation
	}
//1-TIME STEP
//		lograthimcially dsitributes:
	double dt = (1.0/(1.0*alpha_tot))*std::log(1.0/(1.0*ro->get_double()));
	time += dt;

//2-EVENT SLECTION AND ACTION
	double r_int = ro->get_double()*alpha_tot;
	if(r_int <= alpha_in && alpha_in!=0){
		//selecting the type

		int index = search_exceeds_cum<double>(inp,r_int);	
		//adding a monomer
		add_monomer(index);
		return true;
	}		
	else if(r_int <= alpha_out+alpha_in && alpha_out!=0){

		//readjusting the variable
		r_int = r_int - alpha_in;
		//selecting the type
		int index = search_exceeds_cum<double>(outp,r_int);
		//removing a monomer
		remove_monomer(index);
		return true;
	}
	else if(!isZero(alpha_agg)){
	//adjusting the variable
		r_int = (1.0*(r_int -alpha_in-alpha_out))/(R);
	//getting the index and making sure its valid
		if(r_int/cp.get_cum() > 1) throw std::invalid_argument("error with propensities");
		int index = cp.search_exceeds_cum(r_int);
	//getting the two clusters and aggregating
		int row = cp.get_row(index), col = cp.get_col(index);
		if(row >= N_total || col >= N_total || row==col){
			std::cout << alpha_tot << " " << alpha_agg << " "<< cp.get_cum() <<std::endl;
			std::cout << N_total << " " <<row << " " << col << std::endl;

			print_clusters();
			throw std::invalid_argument("out of bounds aggregation");
		}
		aggregate(row,col);
		return true;
	}

	throw std::invalid_argument("shouldnt be here");
	return false;
}
//------------------------------------------------------------------------------
//				OTHER GETTERS
//------------------------------------------------------------------------------
double System::get_time(){
	return time;	
}
int System::get_N_types(){
	return N_types;
}
int System::get_N_total(){
	int sum=0;
	for(int i=0; i<N_types; i++) sum +=N[i];
	return sum;
}
int System::get_N_clusters(){
	return cs.size();
}
int System::get_N(int type){
	return N[type];
}
int System::get_N_m(int type){
	return N_monomers[type];
}
std::vector<std::vector<int>>* System::get_clusters(){
	return &cs;
}

