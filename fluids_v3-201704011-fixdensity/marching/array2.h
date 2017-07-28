#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ARRAY2_H
#define ARRAY2_H

#include "array1.h"
#include <algorithm>
#include <cassert>
#include <vector>

template<class T, class ArrayT=std::vector<T> >
struct Array2
{
   // STL-friendly typedefs

   typedef typename ArrayT::iterator iterator;
   typedef typename ArrayT::const_iterator const_iterator;
   typedef typename ArrayT::size_type size_type;
   typedef long difference_type;
   typedef T& reference;
   typedef const T& const_reference;
   typedef T value_type;
   typedef T* pointer;
   typedef const T* const_pointer;
   typedef typename ArrayT::reverse_iterator reverse_iterator;
   typedef typename ArrayT::const_reverse_iterator const_reverse_iterator;

   // the actual representation

   int ni, nj;
   ArrayT a;

   // the interface

   Array2(void)
      : ni(0), nj(0)
   {}

   Array2(int ni_, int nj_)
      : ni(ni_), nj(nj_), a(ni_*nj_)
   { assert(ni_>=0 && nj>=0); }

   Array2(int ni_, int nj_, ArrayT& a_)
      : ni(ni_), nj(nj_), a(a_)
   { assert(ni_>=0 && nj>=0); }

   Array2(int ni_, int nj_, const T& value)
      : ni(ni_), nj(nj_), a(ni_*nj_, value)
   { assert(ni_>=0 && nj>=0); }

   Array2(int ni_, int nj_, const T& value, size_type max_n_)
      : ni(ni_), nj(nj_), a(ni_*nj_, value, max_n_)
   { assert(ni_>=0 && nj>=0); }

   Array2(int ni_, int nj_, T* data_)
      : ni(ni_), nj(nj_), a(ni_*nj_, data_)
   { assert(ni_>=0 && nj>=0); }

   Array2(int ni_, int nj_, T* data_, size_type max_n_)
      : ni(ni_), nj(nj_), a(ni_*nj_, data_, max_n_)
   { assert(ni_>=0 && nj>=0); }

   template<class OtherArrayT>
   Array2(Array2<T, OtherArrayT>& other)
      : ni(other.ni), nj(other.nj), a(other.a)
   {}

   ~Array2(void)
   {
#ifndef NDEBUG
      ni=nj=0;
#endif
   }

   // single index reference
   const T& operator()(int index) const
   {
	   assert(index>=0 && index < ni*nj);
	   return a[index];
   }

   T& operator()(int index)
   {
	   assert(index>=0 && index < ni*nj);
	   return a[index];
   }

   const T& operator()(int i, int j) const
   {
      assert(i>=0 && i<ni && j>=0 && j<nj);
      return a[i+ni*j];
   }

   T& operator()(int i, int j)
   {
	   if(!(i>=0 && i<ni && j>=0 && j<nj))
		   printf("%d %d\n",i,j);
	  assert(i>=0 && i<ni && j>=0 && j<nj);
      return a[i+ni*j];
   }

   bool operator==(const Array2<T>& x) const
   { return ni==x.ni && nj==x.nj && a==x.a; }

   bool operator!=(const Array2<T>& x) const
   { return ni!=x.ni || nj!=x.nj || a!=x.a; }

   bool operator<(const Array2<T>& x) const
   {
      if(ni<x.ni) return true; else if(ni>x.ni) return false;
      if(nj<x.nj) return true; else if(nj>x.nj) return false;
      return a<x.a;
   }

   bool operator>(const Array2<T>& x) const
   {
      if(ni>x.ni) return true; else if(ni<x.ni) return false;
      if(nj>x.nj) return true; else if(nj<x.nj) return false;
      return a>x.a;
   }

   bool operator<=(const Array2<T>& x) const
   {
      if(ni<x.ni) return true; else if(ni>x.ni) return false;
      if(nj<x.nj) return true; else if(nj>x.nj) return false;
      return a<=x.a;
   }

   bool operator>=(const Array2<T>& x) const
   {
      if(ni>x.ni) return true; else if(ni<x.ni) return false;
      if(nj>x.nj) return true; else if(nj<x.nj) return false;
      return a>=x.a;
   }

   void assign(const T& value)
   { a.assign(value); }

   void assign(int ni_, int nj_, const T& value)
   {
      a.assign(ni_*nj_, value);
      ni=ni_;
      nj=nj_;
   }
    
   void assign(int ni_, int nj_, const T* copydata)
   {
      a.assign(ni_*nj_, copydata);
      ni=ni_;
      nj=nj_;
   }
    
   const T& at(int i, int j) const
   {
      assert(i>=0 && i<ni && j>=0 && j<nj);
      return a[i+ni*j];
   }

   T& at(int i, int j)
   {
      assert(i>=0 && i<ni && j>=0 && j<nj);
      return a[i+ni*j];
   }

   const T& back(void) const
   { 
      assert(a.size());
      return a.back();
   }

   T& back(void)
   {
      assert(a.size());
      return a.back();
   }

   const_iterator begin(void) const
   { return a.begin(); }

   iterator begin(void)
   { return a.begin(); }

   size_type capacity(void) const
   { return a.capacity(); }

   void clear(void)
   {
      a.clear();
      ni=nj=0;
   }

   bool empty(void) const
   { return a.empty(); }

   const_iterator end(void) const
   { return a.end(); }

   iterator end(void)
   { return a.end(); }

   void fill(int ni_, int nj_, const T& value)
   {
      a.fill(ni_*nj_, value);
      ni=ni_;
      nj=nj_;
   }
    
   const T& front(void) const
   {
      assert(a.size());
      return a.front();
   }

   T& front(void)
   {
      assert(a.size());
      return a.front();
   }

   size_type max_size(void) const
   { return a.max_size(); }

   reverse_iterator rbegin(void)
   { return reverse_iterator(end()); }

   const_reverse_iterator rbegin(void) const
   { return const_reverse_iterator(end()); }

   reverse_iterator rend(void)
   { return reverse_iterator(begin()); }

   const_reverse_iterator rend(void) const
   { return const_reverse_iterator(begin()); }

   void reserve(int reserve_ni, int reserve_nj)
   { a.reserve(reserve_ni*reserve_nj); }

   void resize(int ni_, int nj_)
   {
      assert(ni_>=0 && nj_>=0);
      a.resize(ni_*nj_);
      ni=ni_;
      nj=nj_;
   }

   void resize(int ni_, int nj_, const T& value)
   {
      assert(ni_>=0 && nj_>=0);
      a.resize(ni_*nj_, value);
      ni=ni_;
      nj=nj_;
   }

   void set_zero(void)
   { a.set_zero(); }

   size_type size(void) const
   { return a.size(); }

   void swap(Array2<T>& x)
   {
      std::swap(ni, x.ni);
      std::swap(nj, x.nj);
      a.swap(x.a);
   }

   void trim(void)
   { a.trim(); }

   T bilerp(int i, int j, T fx, T fy)
   {  return (1-fx)*((1-fy)*(*this)(i,j)+fy*(*this)(i,j+1))+fx*((1-fy)*(*this)(i+1,j)+fy*(*this)(i+1,j+1)); }

   //! be careful
   //! test it!!
   T monotonic_cubic_interp(int i, int j, T fx, T fy) {

	   if(i<=0||i>=ni-2||j<=0||j>=nj-2)
		   return (1-fx)*((1-fy)*(*this)(i,j)+fy*(*this)(i,j+1))+fx*((1-fy)*(*this)(i+1,j)+fy*(*this)(i+1,j+1));

	   int ineg1 = i-1; if(ineg1<0) ineg1 = 0;
	   int iplus1 = i+1; if(iplus1 >= ni) iplus1 = ni-1;
	   int iplus2 = i+2; if(iplus2 >= ni) iplus2 = ni-1;
	   int jneg1 = j-1; if(jneg1<0) jneg1 = 0;
	   int jplus1 = j+1; if(jplus1 >= nj) jplus1 = nj-1;
	   int jplus2 = j+2; if(jplus2 >= nj) jplus2 = nj-1;

	   T y[4];
	   y[0] = ::monotonic_cubic_interp((*this)(ineg1,jneg1), (*this)(i,jneg1), (*this)(iplus1,jneg1), (*this)(iplus2,jneg1), fx);
	   y[1] = ::monotonic_cubic_interp((*this)(ineg1,j), (*this)(i,j), (*this)(iplus1,j), (*this)(iplus2,j), fx);
	   y[2] = ::monotonic_cubic_interp((*this)(ineg1,jplus1), (*this)(i,jplus1), (*this)(iplus1,jplus1), (*this)(iplus2,jplus1), fx);
	   y[3] = ::monotonic_cubic_interp((*this)(ineg1,jplus2), (*this)(i,jplus2), (*this)(iplus1,jplus2), (*this)(iplus2,jplus2), fx);

	   return ::monotonic_cubic_interp(y[0],y[1],y[2],y[3], fy);
   }

   void copy_to(Array2 &other) const {
	   //std::memcpy(a.a, a, ni*nj*sizeof(T));
	   assert(ni == other.ni && nj == other.nj);
	   for(int i = 0; i < ni*nj; i++)
		   other.a[i] = a[i];
	   //other.a = a;
   }

   T infnorm() const {
	   T r = 0;
	   for(int i = 0; i < ni*nj; i++) {
			if(std::fabs(a[i]) > r)
				r = std::fabs(a[i]);
	   }
	   return r;
   }

};

// some common arrays

typedef Array2<double, Array1<double> >                         Array2d;
typedef Array2<float, Array1<float> >                           Array2f;
typedef Array2<long long, Array1<long long> >                   Array2ll;
typedef Array2<unsigned long long, Array1<unsigned long long> > Array2ull;
typedef Array2<int, Array1<int> >                               Array2i;
typedef Array2<unsigned int, Array1<unsigned int> >             Array2ui;
typedef Array2<short, Array1<short> >                           Array2s;
typedef Array2<unsigned short, Array1<unsigned short> >         Array2us;
typedef Array2<char, Array1<char> >                             Array2c;
typedef Array2<unsigned char, Array1<unsigned char> >           Array2uc;

// and wrapped versions

typedef Array2<double, WrapArray1<double> >                         WrapArray2d;
typedef Array2<float, WrapArray1<float> >                           WrapArray2f;
typedef Array2<long long, WrapArray1<long long> >                   WrapArray2ll;
typedef Array2<unsigned long long, WrapArray1<unsigned long long> > WrapArray2ull;
typedef Array2<int, WrapArray1<int> >                               WrapArray2i;
typedef Array2<unsigned int, WrapArray1<unsigned int> >             WrapArray2ui;
typedef Array2<short, WrapArray1<short> >                           WrapArray2s;
typedef Array2<unsigned short, WrapArray1<unsigned short> >         WrapArray2us;
typedef Array2<char, WrapArray1<char> >                             WrapArray2c;
typedef Array2<unsigned char, WrapArray1<unsigned char> >           WrapArray2uc;

#endif	// ARRAY2_H
