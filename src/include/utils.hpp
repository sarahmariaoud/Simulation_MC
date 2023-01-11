
#ifndef utils_h
#define utils_h

#include <chrono>
#include <string>
#include <sstream>


#define EPS 1.0e-10
inline bool isZero(double d){
    return d<EPS;
}

inline std::string get_time_string(){
time_t rawtime;
  struct tm * timeinfo;
  char buffer [80];

  time (&rawtime);
  timeinfo = localtime (&rawtime);

  strftime (buffer,80,"%F_%H-%M-%S",timeinfo);
  std::string s = buffer;
  return s;
}


template <typename T> int search_exceeds_cum(std::vector<T> arr, T val){
  //cretaing the cumulative sum vector
  std::vector<T> cum_sum;
  cum_sum.resize(arr.size());
  std::partial_sum(arr.begin(),arr.end(),cum_sum.begin());
  //getting the index
  int index = std::lower_bound(cum_sum.begin(),cum_sum.end(),val) - cum_sum.begin();
  if(index>= cum_sum.size()){
   for(int i=0; i<arr.size(); i++){
    std::cout << i << " " << arr[i] <<" " <<cum_sum[i] <<std::endl; 
   }
   std::cout << " VAL " << val << " INDEX " << index << std::endl; 
   throw std::invalid_argument("the value exceed the matrix");
  }

  return  (index);

  
}

inline std::string tostr(double value) {
    std::string out = std::string("");
    if(value < 0 ){
      value = -1*value;
      out = std::string("n");
    }
    std::ostringstream s;
    s.precision(std::numeric_limits<double>::digits10);
    s << value;
    out = out+s.str();
    return out;
}


#endif