///with the help of this guy:
//https://stackoverflow.com/questions/2808398/easily-measure-elapsed-time


#ifndef tick_tock_h
#define tick_tock_h

#include <chrono>
#include <thread>
#include <assert.h>
using namespace std::chrono_literals;


template <class DT = std::chrono::milliseconds, class ClockT = std::chrono::steady_clock>
class Timer
{
    using timep_t = decltype(ClockT::now());
    
    timep_t _start = ClockT::now();
    timep_t _end = {};

public:
    void tick() { 
        _end = timep_t{};
        _start = ClockT::now(); 
    }
    
    void tock() {
        _end = ClockT::now(); 
    }
    
    template <class duration_t = DT> auto duration() const { 
        // Use gsl_Expects if your project supports it.
        assert(_end != timep_t{} && "Timer must toc before reading the time"); 
        return std::chrono::duration_cast<duration_t>(_end - _start).count(); 
    }
};

#endif