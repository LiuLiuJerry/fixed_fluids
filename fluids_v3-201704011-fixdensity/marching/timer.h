#if defined(_MSC_VER)
#pragma once
#endif

#ifndef TIMER_H
#define TIMER_H

#include <time.h>

// Class used to time the duration of events.
// This class can be used as a stopwatch.
class Timer {
public:
   // Construct a new timer, initialized to zero.
   Timer(); 
   // Reset the timer to zero.
   void Reset();
   // Clear the timers memory.
   void Clear();
   // Start the timer
   void Start();
   // Stop the timer
   void Stop();
   // Returns the current value of the timer, expressed in msecs.
   double Value();
   // Returns the current value of the timer, expressed in secs.
   double GetSeconds();
   // Returns the current value of the timer, expressed in Hz.
   double FPS();
   double AverageFPS();
   // Returns the average time over all runs of the timer from the last reset.
   double Average();
protected:
   clock_t begin;	// Begin clock tick
   clock_t end;	// End clock tick.
private:
   double total;	// Total time over all runs
   double nb;		// Number of runs since last reset
   double msec;	// Time in msec
};


#endif	// TIMER_H
