#ifndef _MMATH_H_
#define _MMATH_H_

#include "kernel.h"
#include "twister.h"
# define Pi 3.14

extern TWISTER Rng;

inline double clamp (double val, double Lo, double Hi)
  {
  if (val < Lo) return Lo;
  if (val > Hi) return Hi;
  return val;
  }

inline double cycloidal (double x)
  {
  return sin (x*2*Pi);
  }

inline double trianglewave (double x)
  {
  double offset = fmod(x, double(1));
  if (offset < 0)
    offset += 1;
  if (offset > 0.5)
    offset = 1 - offset;
  return offset + offset;
  }

inline double lerp (double t, double lo, double hi)
  {
  return lo + t*(hi-lo);
  }

inline double scurve (double val)
  {
  return (val*val*(3. - 2.*val));  
  }

inline double posmod(double x, double y)
  {
  double ret = fmod(x, y);
  if (ret < 0.0) ret += y;
  return ret / y;
  }

inline double dblrand (double min, double max)
  {
  return (max-min) * (Rng.f8rand()) + min;
  }

#endif
