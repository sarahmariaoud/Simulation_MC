

#ifndef random_object_h
#define random_object_h

#include <random>
#include <chrono>

class RandomObject{
    std::default_random_engine generator;
   
public:
    RandomObject();
    RandomObject(int seed);
    double get_double();
    int get_int(int start, int end);
};


inline RandomObject::RandomObject(){
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());    
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());    
}

inline RandomObject::RandomObject(int seed){

    generator.seed(seed);    
}
inline double RandomObject::get_double(){
    std::uniform_real_distribution<double> distribution;
    return distribution(generator);
}
inline int RandomObject::get_int(int start, int end){
    std::uniform_int_distribution<int> distribution(start, end);
    return distribution(generator);

}




#endif