#if defined(_MSC_VER)
#pragma once
#endif

#ifndef UTIL_H
#define UTIL_H

// C_plus_plus header
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <cmath>
#include <string.h>
#include <string>
#include <assert.h>
#include <algorithm>
#include <malloc.h>
#include <float.h>
#include <time.h>
#include <map>


//#include <gl\glut.h>

#ifdef WIN32

// I hate warning...
#pragma warning (disable: 4996)
#pragma warning (disable: 4018)

#endif


#ifdef WIN32
#undef min
#undef max
#endif

using std::vector;
using std::set;
using std::string;
using std::max;
using std::min;
using std::sort;
using std::swap;

// Console print function
#include "console.h"

#define isnan _isnan
/*
#define isnan(x) ((x) != (x))
*/
#define isinf(f) (!_finite((f)))

//#define uint unsigned int 
typedef unsigned int		uint;
#define int8_t __int8
#define uint8_t unsigned __int8
#define int16_t __int16
#define uint16_t unsigned __int16
#define int32_t __int32
#define uint32_t unsigned __int32
#define int64_t __int64
#define uint64_t unsigned __int64

// Safe Delete Macro, from DXUT.h
#ifndef SAFE_DELETE
#define SAFE_DELETE(p)       { if (p) { delete (p);     (p)=NULL; } }
#endif
#ifndef SAFE_DELETE_ARRAY
#define SAFE_DELETE_ARRAY(p) { if (p) { delete[] (p);   (p)=NULL; } }
#endif

// use custom random number generator
#define USE_CUSTOM_RNG

#ifdef NDEBUG
#define Assert(expr) ((void)0)
#else
#define Assert(expr) \
	((expr) ? (void)0 : \
	Fatal("Assertion \"%s\" failed in %s, line %d", \
#expr, __FILE__, __LINE__))
#endif

// Global Constants
#define FS_VERSION "0.0.4"
#ifdef M_PI
#undef M_PI
#endif
#define M_PI		3.14159265358979323846
#define INV_PI		0.31830988618379067154
#define INV_TWOPI	0.15915494309189533577
#define INV_FOURPI	0.07957747154594766788

static const float OneMinusEpsilon=0.9999999403953552f;

#ifndef INFINITY
#define INFINITY FLT_MAX
#endif

// memory allocation

// temporary memory allocation, automatically delete, used as local memory
#define ALLOCA(TYPE, COUNT) (TYPE*)alloca((COUNT)*sizeof*TYPE))

//! When use template function, specify the type explicitly!!!
// especially min and max...sometimes it doesn't work. e.g. interpPhi function in fs_pls.h

// util functions
template<class T>
inline T sqr(const T& x)
{ return x*x; }

template<class T>
inline T cube(const T& x)
{ return x*x*x; }

template<class T>
inline T min(T a1, T a2, T a3)
{ return min(a1, min(a2, a3)); }

template<class T>
inline T min(T a1, T a2, T a3, T a4)
{ return min(min(a1, a2), min(a3, a4)); }

template<class T>
inline T min(T a1, T a2, T a3, T a4, T a5)
{ return min(min(a1, a2), min(a3, a4), a5); }

template<class T>
inline T min(T a1, T a2, T a3, T a4, T a5, T a6)
{ return min(min(a1, a2), min(a3, a4), min(a5, a6)); }

template<class T>
inline T max(T a1, T a2, T a3)
{ return max(a1, max(a2, a3)); }

template<class T>
inline T max(T a1, T a2, T a3, T a4)
{ return max(max(a1, a2), max(a3, a4)); }

template<class T>
inline T max(T a1, T a2, T a3, T a4, T a5)
{ return max(max(a1, a2), max(a3, a4),  a5); }

template<class T>
inline T max(T a1, T a2, T a3, T a4, T a5, T a6)
{ return max(max(a1, a2), max(a3, a4),  max(a5, a6)); }

template<class T>
inline void minmax(T a1, T a2, T& amin, T& amax)
{
	if(a1<a2){
		amin=a1;
		amax=a2;
	}else{
		amin=a2;
		amax=a1;
	}
}

template<class T>
inline void minmax(T a1, T a2, T a3, T& amin, T& amax)
{
	if(a1<a2){
		if(a1<a3){
			amin=a1;
			if(a2<a3) amax=a3;
			else amax=a2;
		}else{
			amin=a3;
			if(a1<a2) amax=a2;
			else amax=a1;
		}
	}else{
		if(a2<a3){
			amin=a2;
			if(a1<a3) amax=a3;
			else amax=a1;
		}else{
			amin=a3;
			amax=a1;
		}
	}
}

template<class T>
inline void minmax(T a1, T a2, T a3, T a4, T& amin, T& amax)
{
	if(a1<a2){
		if(a3<a4){
			amin=min(a1,a3);
			amax=max(a2,a4);
		}else{
			amin=min(a1,a4);
			amax=max(a2,a3);
		}
	}else{
		if(a3<a4){
			amin=min(a2,a3);
			amax=max(a1,a4);
		}else{
			amin=min(a2,a4);
			amax=max(a1,a3);
		}
	}
}

template<class T>
inline void minmax(T a1, T a2, T a3, T a4, T a5, T& amin, T& amax)
{
	//@@@ the logic could be shortcircuited a lot!
	amin=min(a1,a2,a3,a4,a5);
	amax=max(a1,a2,a3,a4,a5);
}

template<class T>
inline void minmax(T a1, T a2, T a3, T a4, T a5, T a6, T& amin, T& amax)
{
	//@@@ the logic could be shortcircuited a lot!
	amin=min(a1,a2,a3,a4,a5,a6);
	amax=max(a1,a2,a3,a4,a5,a6);
}

template<class T>
inline void update_minmax(T a1, T& amin, T& amax)
{
	if(a1<amin) amin=a1;
	else if(a1>amax) amax=a1;
}

template<class T>
inline void sort(T &a, T &b, T &c)
{
	T temp;
	if(a<b){
		if(a<c){
			if(c<b){ // a<c<b
				temp=c;c=b;b=temp;
			} // else: a<b<c
		}else{ // c<a<b
			temp=c;c=b;b=a;a=temp;
		}
	}else{
		if(b<c){
			if(a<c){ //b<a<c
				temp=b;b=a;a=temp;
			}else{ // b<c<a
				temp=b;b=c;c=a;a=temp;
			}
		}else{ // c<b<a
			temp=c;c=a;a=temp;
		}
	}
}

template<class T>
inline T clamp(T a, T lower, T upper)
{
	if(a<lower) return lower;
	else if(a>upper) return upper;
	else return a;
}

// only makes sense with T=float or double
template<class T>
inline T smooth_step(T r)
{
	if(r<0) return 0;
	else if(r>1) return 1;
	return r*r*r*(10+r*(-15+r*6));
}

// only makes sense with T=float or double
template<class T>
inline T smooth_step(T r, T r_lower, T r_upper, T value_lower, T value_upper)
{ return value_lower + smooth_step((r-r_lower)/(r_upper-r_lower)) * (value_upper-value_lower); }

// only makes sense with T=float or double
template<class T>
inline T ramp(T r)
{ return smooth_step((r+1)/2)*2-1; }

#ifdef WIN32
inline int lround(double x)
{
	if(x>0)
		return (x-floor(x)<0.5) ? (int)floor(x) : (int)ceil(x);
	else
		return (x-floor(x)<=0.5) ? (int)floor(x) : (int)ceil(x);
}

inline double remainder(double x, double y)
{
	return x-std::floor(x/y+0.5)*y;
}
#endif

inline unsigned int round_up_to_power_of_two(unsigned int n)
{
	int exponent=0;
	--n;
	while(n){
		++exponent;
		n>>=1;
	}
	return 1<<exponent;
}

inline unsigned int round_down_to_power_of_two(unsigned int n)
{
	int exponent=0;
	while(n>1){
		++exponent;
		n>>=1;
	}
	return 1<<exponent;
}

// Transforms even the sequence 0,1,2,3,... into reasonably good random numbers 
// Challenge: improve on this in speed and "randomness"!
// This seems to pass several statistical tests, and is a bijective map (of 32-bit unsigned ints)
inline unsigned int randhash(unsigned int seed)
{
	unsigned int i=(seed^0xA3C59AC3u)*2654435769u;
	i^=(i>>16);
	i*=2654435769u;
	i^=(i>>16);
	i*=2654435769u;
	return i;
}

// the inverse of randhash
inline unsigned int unhash(unsigned int h)
{
	h*=340573321u;
	h^=(h>>16);
	h*=340573321u;
	h^=(h>>16);
	h*=340573321u;
	h^=0xA3C59AC3u;
	return h;
}

// returns repeatable stateless pseudo-random number in [0,1]
inline double randhashd(unsigned int seed)
{ return randhash(seed)/(double)UINT_MAX; }
inline float randhashf(unsigned int seed)
{ return randhash(seed)/(float)UINT_MAX; }

// returns repeatable stateless pseudo-random number in [a,b]
inline double randhashd(unsigned int seed, double a, double b)
{ return (b-a)*randhash(seed)/(double)UINT_MAX + a; }
inline float randhashf(unsigned int seed, float a, float b)
{ return ( (b-a)*randhash(seed)/(float)UINT_MAX + a); }

inline int intlog2(int x)
{
	int exp=-1;
	while(x){
		x>>=1;
		++exp;
	}
	return exp;
}

template<class T>
inline void get_barycentric(T x, int& i, T& f, int i_low, int i_high)
{
	T s=std::floor(x);
	i=(int)s;
	if(i<i_low){
		i=i_low;
		f=0;
	}else if(i>i_high-2){
		i=i_high-2;
		f=1;
	}else
		f=(T)(x-s);
}

template<class S, class T>
inline S lerp(const S& value0, const S& value1, T f)
{ return (1-f)*value0 + f*value1; }

template<class S, class T>
inline S bilerp(const S& v00, const S& v10, 
	const S& v01, const S& v11, 
	T fx, T fy)
{ 
	return lerp(lerp(v00, v10, fx),
		lerp(v01, v11, fx), 
		fy);
}

template<class S, class T>
inline S trilerp(const S& v000, const S& v100,
	const S& v010, const S& v110,
	const S& v001, const S& v101,  
	const S& v011, const S& v111,
	T fx, T fy, T fz) 
{
	return lerp(bilerp(v000, v100, v010, v110, fx, fy),
		bilerp(v001, v101, v011, v111, fx, fy),
		fz);
}

template<class S, class T>
inline S quadlerp(const S& v0000, const S& v1000,
	const S& v0100, const S& v1100,
	const S& v0010, const S& v1010,  
	const S& v0110, const S& v1110,
	const S& v0001, const S& v1001,
	const S& v0101, const S& v1101,
	const S& v0011, const S& v1011,  
	const S& v0111, const S& v1111,
	T fx, T fy, T fz, T ft) 
{
	return lerp(trilerp(v0000, v1000, v0100, v1100, v0010, v1010, v0110, v1110, fx, fy, fz),
		trilerp(v0001, v1001, v0101, v1101, v0011, v1011, v0111, v1111, fx, fy, fz),
		ft);
}

// f should be between 0 and 1, with f=0.5 corresponding to balanced weighting between w0 and w2
template<class T>
inline void quadratic_bspline_weights(T f, T& w0, T& w1, T& w2)
{
	w0=T(0.5)*sqr(f-1);
	w1=T(0.75)-sqr(f-T(0.5));;
	w2=T(0.5)*sqr(f);
}

// f should be between 0 and 1
template<class T>
inline void cubic_interp_weights(T f, T& wneg1, T& w0, T& w1, T& w2)
{
	T f2(f*f), f3(f2*f);
	wneg1=-T(1./3)*f+T(1./2)*f2-T(1./6)*f3;
	w0=1-f2+T(1./2)*(f3-f);
	w1=f+T(1./2)*(f2-f3);
	w2=T(1./6)*(f3-f);
}

template<class S, class T>
inline S cubic_interp(const S& value_neg1, const S& value0, const S& value1, const S& value2, T f)
{
	T wneg1, w0, w1, w2;
	cubic_interp_weights(f, wneg1, w0, w1, w2);
	return wneg1*value_neg1 + w0*value0 + w1*value1 + w2*value2;
}

// monotonic cubic interpolation
// "Visual Simulation of Smoke" by Fedkiw
template<class S, class T>
inline S monotonic_cubic_interp(const S& value_neg1, const S& value0, const S& value1, const S& value2, T f)
{
	/*S delta_k = (value1 - value0);
	S d_k = (value1 - value_neg1) / 2;
	S d_k1 = (value2 - value0) / 2;

	if(delta_k == 0) {
		d_k = d_k1 = 0;
	}
	else {
		if(sign(d_k) != sign(delta_k))
			d_k = 0;
		if(sign(d_k1) != sign(delta_k))
			d_k1 = 0;
	}

	S a0 = value0;
	S a1 = d_k;
	S a2 = 3*delta_k - 2*d_k - d_k1;
	S a3 = d_k + d_k1 - 2*delta_k;

	return f*(a1+f*(a2+a3*f))+a0;*/

	// Bridson' Fluid Simulation Page 38
	T f2 = f*f, f3 = f*f2;

	S ret= value_neg1*(-0.5*f+f2-0.5*f3)+value0*(1-2.5*f2+1.5*f3)+value1*(0.5*f+2*f2-1.5*f3)+ value2*(-0.5*f2+0.5*f3);
	ret = clamp(ret, min<S>(value0,value1),max<S>(value0,value1));
	return ret;
}

template<class T>
void zero(std::vector<T>& v)
{ for(int i=(int)v.size()-1; i>=0; --i) v[i]=0; }

template<class T>
T abs_max(const std::vector<T>& v)
{
	T m=0;
	for(int i=(int)v.size()-1; i>=0; --i){
		if(std::fabs(v[i])>m)
			m=std::fabs(v[i]);
	}
	return m;
}

template<class T>
bool contains(const std::vector<T>& a, T e)
{
	for(unsigned int i=0; i<a.size(); ++i)
		if(a[i]==e) return true;
	return false;
}

template<class T>
void add_unique(std::vector<T>& a, T e)
{
	for(unsigned int i=0; i<a.size(); ++i)
		if(a[i]==e) return;
	a.push_back(e);
}

template<class T>
void insert(std::vector<T>& a, unsigned int index, T e)
{
	a.push_back(a.back());
	for(unsigned int i=(unsigned int)a.size()-1; i>index; --i)
		a[i]=a[i-1];
	a[index]=e;
}

template<class T>
void erase(std::vector<T>& a, unsigned int index)
{
	for(unsigned int i=index; i<a.size()-1; ++i)
		a[i]=a[i+1];
	a.pop_back();
}

template<class T>
void erase_swap(std::vector<T>& a, unsigned int index)
{
	for(unsigned int i=index; i<a.size()-1; ++i)
		swap(a[i], a[i+1]);
	a.pop_back();
}

template<class T>
void erase_unordered(std::vector<T>& a, unsigned int index)
{
	a[index]=a.back();
	a.pop_back();
}

template<class T>
void erase_unordered_swap(std::vector<T>& a, unsigned int index)
{
	swap(a[index], a.back());
	a.pop_back();
}

template<class T>
void find_and_erase_unordered(std::vector<T>& a, const T& doomed_element)
{
	for(unsigned int i=0; i<a.size(); ++i)
		if(a[i]==doomed_element){
			erase_unordered(a, i);
			return;
		}
}

template<class T>
void replace_once(std::vector<T>& a, const T& old_element, const T& new_element)
{
	for(unsigned int i=0; i<a.size(); ++i)
		if(a[i]==old_element){
			a[i]=new_element;
			return;
		}
}

template<class T>
inline int sign(T x) {
	if (x>(T)0) return 1;
	else return -1;
}

template<class T>
inline T absmin2 (T x, T y) {
	if (fabs(x) < fabs(y)) return x;
	return y;
}

template<class T>
inline T absmax2 (T x, T y) {
	if (fabs(x) > fabs(y)) return x;
	return y;
}

template <class T> 
inline void Matrix_Inverse_3(T *A, T *R) { //R=inv(A)
	R[0]=A[4]*A[8]-A[7]*A[5];
	R[1]=A[7]*A[2]-A[1]*A[8];
	R[2]=A[1]*A[5]-A[4]*A[2];
	R[3]=A[5]*A[6]-A[3]*A[8];
	R[4]=A[0]*A[8]-A[2]*A[6];
	R[5]=A[2]*A[3]-A[0]*A[5];
	R[6]=A[3]*A[7]-A[4]*A[6];
	R[7]=A[1]*A[6]-A[0]*A[7];
	R[8]=A[0]*A[4]-A[1]*A[3];
	T inv_det=(T)1/(A[0]*R[0]+A[3]*R[1]+A[6]*R[2]);
	if(isnan(inv_det))
		Fatal("Error in matrix inverse\n");
	for(int i=0; i<9; i++)
		R[i]*=inv_det;
}

template <class T> 
inline void Matrix_Vector_Product_3(T *A, T *x, T *r) {	//r=A*x
	r[0]=A[0]*x[0]+A[1]*x[1]+A[2]*x[2];
	r[1]=A[3]*x[0]+A[4]*x[1]+A[5]*x[2];
	r[2]=A[6]*x[0]+A[7]*x[1]+A[8]*x[2];
}

template <class T> 
inline void Matrix_Transpose_3(T *A, T *R) //R=A'
{
	memcpy(R, A, sizeof(T)*9);
	swap(R[1], R[3]);
	swap(R[2], R[6]);
	swap(R[5], R[7]);
}

template <class T>
inline T MagnitudeArray3(T *x)
{
	return sqrtf(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

template <class T>
inline void Matrix_Product_3(T *A, T *B, T *R)	{	//R=A*B
	R[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];
	R[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];
	R[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];
	R[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];
	R[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];
	R[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];
	R[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];
	R[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];
	R[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
}

template <class T> 
inline void Cross(T *x, T *y, T *r) {
	r[0]= x[1]*y[2]-y[1]*x[2];
	r[1]=-x[0]*y[2]+y[0]*x[2];
	r[2]= x[0]*y[1]-y[0]*x[1];
}

template <class T> inline T Normalize(T *x, T *r=0) {
	if(r==0)	r=x;
	T m=Magnitude(x);
	if(m<1e-14f)	{printf("ERROR: vector cannot be normalized.\n"); return m;}
	T inv_m=1/m;
	r[0]=x[0]*inv_m;
	r[1]=x[1]*inv_m;
	r[2]=x[2]*inv_m;
	return m;
}

template <class T>
inline T Magnitude_Squared(T *x) {
	return x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
}

template <class T>
inline T Magnitude(T *x) {
	return sqrtf(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

template <class T> 
inline T Dot(T *x, T *y, int number=3){
	T ret=0;
	for(int i=0; i<number; i++)	ret+=x[i]*y[i];
	return ret;
}

#endif	// UTIL_H