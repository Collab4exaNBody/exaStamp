/// @file
/// @brief

#ifndef __CHRONO_HPP_INCLUDED
#define __CHRONO_HPP_INCLUDED


#ifdef __use_cycle_h


extern "C"
{
#include "cycle.h"
#include <unistd.h>
}


#include <sys/time.h>


/// @brief Class to measure times
class Chrono {

public:

  Chrono()
  {
    this->got_first = 0;
	this->got_last = 0;  
	this->ticks_per_sec = 0.0;  
  }
  
  ~Chrono() {}

  void tic() { 
	if( !got_first )
	{
		got_first = 1;
		gettimeofday(&firstcall, NULL);
	}
    start = getticks(); 
    
    //fprintf(stderr, "TIC %llu\n", start );
  }

  void toc() {
	//We use the first delta to init the
	//ticks_per_sec value
    if( !got_last )
	{
		got_last = 1;
		
		ticks now = getticks();
		
		if( (now - start) < 1e6 )
		{
			fprintf(stderr, "Calibrating timestamp ... \n");
			sleep( 1 );
			now = getticks();
		}
		
		gettimeofday(&lastcall, NULL);
		
		/* Now compute tps */
		double start_time = (firstcall.tv_usec) * 1.0e-06 + (firstcall.tv_sec) * 1.0;
		double end_time = (lastcall.tv_usec) * 1.0e-06 + (lastcall.tv_sec) * 1.0;
		double tps = ( (double)(now - start) ) / ( end_time - start_time );
		
		ticks_per_sec = tps;
	}
    
    
    stop = getticks();
    //fprintf(stderr, "TOC %llu D %llu\n", stop, stop - start );
  }

  double elapsed() {
	
	//fprintf(stderr, "\n=> %g / %g\n" , ((double)(stop - start))/ticks_per_sec , ticks_per_sec);
	
    return ((double)(stop - start))/ticks_per_sec;
  }

private:  

  ticks start, stop;
  timeval firstcall, lastcall;
  int got_first, got_last;
  double ticks_per_sec;
};


#elif  defined(__use_c_chrono)

#include <chrono>

/// @brief Class to measure
class Chrono {

  typedef std::chrono::steady_clock Clock_t;

public:

  Chrono() {}
  ~Chrono() {}

  void tic() { 
    start = Clock_t::now(); 
  }

  void toc() {
    stop  = Clock_t::now();
  }

  double elapsed() { 
    return std::chrono::duration_cast< std::chrono::microseconds >(stop-start).count() * 1.0e-06;
  }

private:

  Clock_t::time_point start, stop;

};


#else


#include <sys/time.h>

/// @brief Class to measure times
class Chrono {

public:

  /// @brief Constructor
  Chrono() {}
  
  /// @brief Destructor
  ~Chrono() {}

  /// @brief Get a first time
  void tic() { 
    gettimeofday(&start, NULL); 
  }

  /// @brief Get the second measure of time
  void toc() { 
    gettimeofday(&stop , NULL);
  }

  /// @brief Compute time between tic and toc
  /// @return Value of the chrono
  double elapsed() { 
    return (stop.tv_usec - start.tv_usec) * 1.0e-06 + (stop.tv_sec - start.tv_sec) * 1.0;
  }

private:  

  timeval start; ///< first time (tic)
  timeval stop;  ///< second time (toc)

};



#endif




#endif //  __CHRONO_HPP_INCLUDED
