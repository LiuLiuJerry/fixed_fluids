#if defined(_MSC_VER)
#pragma once
#endif

#ifndef TYPES_H
#define TYPES_H

#include "array2.h"
#include "array3.h"
#include "vec.h"

/************************************************************************/
/*					Type declaration                                    */
/************************************************************************/

typedef float								Real;
typedef Real								real;
typedef Vec<2,Real>							Vec2r;
typedef Array2<Real, Array1<Real>>			Array2r;

typedef Vec<3,Real>							Vec3r;
typedef Array3<Real, Array1<Real>>			Array3r;

typedef Vec<4,Real>							Vec4r;

typedef unsigned char uchar;
typedef signed char schar;

#endif	// TYPES_H