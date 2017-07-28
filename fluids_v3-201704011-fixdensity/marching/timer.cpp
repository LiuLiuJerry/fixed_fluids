#include "util.h"
#include "timer.h"

Timer::Timer() {
   total = 0;
   nb = 0;
   Reset();
}

void Timer::Clear() {
   total = 0;
   nb = 0;
   Reset();
}

void Timer::Reset() {
   msec = 0;
   end = begin;
}

void Timer::Start() {
   begin = clock();
}

void Timer::Stop() {
   end = clock();
   double time = end - begin;
   msec += time;
   begin = end;
}

double Timer::Value() {
   double result = msec;
   total += result;
   nb+=1;
   return result;
}

double Timer::GetSeconds() {
	return Value() / 1000.0;
}

double Timer::FPS() {
   return 1000.0 / Value();
}

double Timer::AverageFPS() {
   return 1000.0 / Average();
}

double Timer::Average() {
   return total / nb;
}
