////////////////////////////////////////////////////////////////////////////////////////
//					LOWER TRIANGLE STRUCTURE
////////////////////////////////////////////////////////////////////////////////////////
//
//	https://math.stackexchange.com/questions/646117/how-to-find-a-function-mapping-matrix-indices
//
//	possible modification:
//		
//		- create overflow protection
//		- instead of creating maps should the indices be calculated each time
//			what is more efficient (one memory and one computation)
//			can test later
////////////////////////////////////////////////////////////////////////////////////////



#ifndef lower_triangle_h
#define lower_triangle_h

#include <iostream>
#include <map>
#include <math.h>
#include <numeric>
#include <vector>


////////////////////////////////////////////////////////////////////////////////////////
//						CLASS DEFINITION
////////////////////////////////////////////////////////////////////////////////////////


template <typename T> class LowerTriangle{
	
public:
	int dim; 	//dimension 
	int size;	//total size
	std::vector<T> arr;
	
	//have the total value of the lower triangle already
	T cumulative; //this is the cumulative of lower triangle including diagonal only


	//to create an easy search table
	std::map<int, std::pair<int,int> > map_index;
	std::map<std::pair<int,int>,int > map_row_col;

public:
	LowerTriangle(int dim_);
	LowerTriangle(const LowerTriangle & lt);
	~LowerTriangle();
	//main setters and getters
	void set(int r, int c, T val);
	void set(int i, T val);
	T get(int r, int c);
	T get(int i);
	//changing the size
	void remove(int n);
	void add();

	//retrieving some variables
	int get_dim();
	int get_size();
	int get_index(int row, int col);
	int get_row(int index);
	int get_col(int index);
	T get_cum();
	//function to get index of first element that exceeds the cumulative sum
	int search_exceeds_cum(T value);

	//printing functions
	void print_array();
	void print_lower_triangle();
	void print_indexed();


private:
	int get_index_from_row_col(int row, int col);
	int get_row_from_index(int index);
	int get_col_from_index(int index);

};

////////////////////////////////////////////////////////////////////////////////////////
//						Constructor/Deconstructor
////////////////////////////////////////////////////////////////////////////////////////
template <typename T> LowerTriangle<T>::LowerTriangle(int dim_): dim(dim_){
	//determining size
	size= (int) ((dim*(dim+1))*0.5);
	//creating the matrix and search table
	arr.resize(size,0);
	for(int i=0; i<size ; i++){
		 map_index[i] = {get_row_from_index(i),get_col_from_index(i)};
		 map_row_col[{get_row_from_index(i),get_col_from_index(i)}] = i;
	}
	cumulative =0;
}
template<typename T> LowerTriangle<T>::LowerTriangle(const LowerTriangle  & lt):
dim(lt.dim),size(lt.size), arr(lt.arr),map_index(lt.map_index),map_row_col(lt.map_row_col){}
template <typename T> LowerTriangle<T>::~LowerTriangle(){}
////////////////////////////////////////////////////////////////////////////////////////
//						Constructor/Deconstructor
////////////////////////////////////////////////////////////////////////////////////////

template <typename T>  void LowerTriangle<T>::remove(int n){
	
	if(n>=dim) throw std::invalid_argument("cant remove that element matrix size exceeded");

	for(int i=0; i<dim; i++){
		cumulative -= arr[get_index(n,i)-i];
		arr.erase(arr.begin()+ (get_index(n,i)-i));
	}	
	size = size - dim; 
	dim = dim-1;

	//No need to resize the index maps we can just keep them and increse them when needed
	
}
template <typename T> void LowerTriangle<T>::add(){
	//getting the new dimensions
	int dim_new = dim+1;
	int size_new = size+ dim_new;

	arr.resize(size_new,0);

	//check if we need to increase the maps	
	if(size_new>size){
		for(int i=size; i<size_new; i++){
			map_index[i] = {get_row_from_index(i),get_col_from_index(i)};
		 	map_row_col[{get_row_from_index(i),get_col_from_index(i)}] = i;
		}
	}
	dim = dim_new;
	size = size_new;

}

////////////////////////////////////////////////////////////////////////////////////////
//						getters and setters of array
////////////////////////////////////////////////////////////////////////////////////////

template <typename T> T LowerTriangle<T>::get(int r, int c){
	if(r>= dim || c>= dim) throw std::invalid_argument("exceeds dim");
	return arr[map_row_col[{std::max(r,c),std::min(r,c)}]];	
}
template <typename T> T LowerTriangle<T>::get(int i){
	if(i>=size) throw std::invalid_argument("exceeds size");
	return arr[i];	
}
template <typename T> void LowerTriangle<T>::set(int r, int c, T val){
	if(r>= dim || c>= dim) throw std::invalid_argument("exceeds dim");
	//editing cumulative
	cumulative -= arr[map_row_col[{std::max(r,c),std::min(r,c)}]];
	cumulative += val;
	//setting the array
	arr[map_row_col[{std::max(r,c),std::min(r,c)}]] = val;	
}
template <typename T> void LowerTriangle<T>::set(int i ,T val){
	if(i>=size) throw std::invalid_argument("exceeds size");
	//editing the cumulative
	cumulative -= arr[i];
	cumulative += val;
	arr[i] = val;	
}


////////////////////////////////////////////////////////////////////////////////////////
//						getters and setters of object
////////////////////////////////////////////////////////////////////////////////////////
template <typename T> int LowerTriangle<T>::get_dim() {return dim;}
template <typename T> int LowerTriangle<T>::get_size() {return size;}
template <typename T> int LowerTriangle<T>::get_index(int r, int c){ return map_row_col[{std::max(r,c),std::min(r,c)}]; }
template <typename T> int LowerTriangle<T>::get_row(int index) {return map_index[index].first;}
template <typename T> int LowerTriangle<T>::get_col(int index) {return map_index[index].second;}
template <typename T> T LowerTriangle<T>::get_cum(){ return cumulative;}

////////////////////////////////////////////////////////////////////////////////////////
//						searching algorithm
////////////////////////////////////////////////////////////////////////////////////////
template <typename T> int  LowerTriangle<T>::search_exceeds_cum(T val){
	//creating the cumulative sum matrix
	std::vector<T> cum_sum;
	cum_sum.resize(arr.size());
    std::partial_sum(arr.begin(),arr.end(),cum_sum.begin());
    //getting the index
    int index = std::lower_bound(cum_sum.begin(), cum_sum.end(),val)- cum_sum.begin();

    if(index>= cum_sum.size()) throw std::invalid_argument("the value exceed the matrix");
    return index;
   
}
////////////////////////////////////////////////////////////////////////////////////////
//						printing the array
////////////////////////////////////////////////////////////////////////////////////////

template <typename T> void LowerTriangle<T>::print_array(){
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			std::cout << get(i,j) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
template <typename T> void LowerTriangle<T>::print_lower_triangle(){
	for(int i=0; i<dim; i++){
		for(int j=0; j<=i; j++){
			std::cout << get(i,j) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
template <typename T> void LowerTriangle<T>::print_indexed(){
	T cum_sum =0;
	for(int i=0; i<size; i++){
		cum_sum += arr[i];
		std::cout << i  << " ( " << map_index[i].first << " , " << map_index[i].second << ") --> " << arr[i] <<"\t" << cum_sum <<std::endl;

	}
}
////////////////////////////////////////////////////////////////////////////////////////
//						private functions for initializing the matrix
////////////////////////////////////////////////////////////////////////////////////////

template <typename T> int LowerTriangle<T>::get_index_from_row_col(int row, int col){
	return ((0.5*std::max(row,col)*(std::max(row,col)+1))+(std::min(row,col)+1))-1;
}

template <typename T> int LowerTriangle<T>::get_row_from_index(int index){
	return (int) (std::ceil(sqrt(2 *(index+1) + 0.25) - 0.5)-1);
}

template <typename T> int LowerTriangle<T>::get_col_from_index(int index){
	int row = get_row_from_index(index);
	return ((index+1) - ((row*(row+1))/2) )-1; 
}

////////////////////////////////////////////////////////////////////////////////////////
//						END OF HEADER FILE
////////////////////////////////////////////////////////////////////////////////////////

#endif

