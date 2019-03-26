/**
 * \file OpenMPTimer.h
 * \brief A simple timer class.
 *
 * \author Pierre Kestener
 * \date 29 Oct 2010
 *
 */
#ifndef OPENMP_TIMER_H_
#define OPENMP_TIMER_H_

#include <omp.h>

/**
 * \brief a simple Timer class.
 * If MPI is enabled, should we use MPI_WTime instead of gettimeofday (?!?)
 */
class OpenMPTimer
{
public:
  /** default constructor, timing starts rightaway */
  OpenMPTimer() {
    start_time = 0.0;
    total_time = 0.0;
    start();
  };
    
  OpenMPTimer(double t) {
      start_time = 0;
      total_time = t;
  };

  OpenMPTimer(OpenMPTimer const& aTimer)  : start_time(aTimer.start_time), total_time(aTimer.total_time) {};
  virtual ~OpenMPTimer();

  /** start time measure */
  virtual void start() { start_time = omp_get_wtime(); };
  
  /** stop time measure and add result to total_time */
  virtual void stop() {
    double now = omp_get_wtime();
    total_time += (now-start_time);
  };

  /** return elapsed time in seconds (as stored in total_time) */
  virtual double elapsed() const {  return total_time; };

protected:
  double    start_time;

  /** store total accumulated timings */
  double    total_time;

}; // class OpenMPTimer


#endif // OPENMP_TIMER_H_
