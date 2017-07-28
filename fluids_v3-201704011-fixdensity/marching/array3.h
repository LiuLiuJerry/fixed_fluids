#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ARRAY3_H
#define ARRAY3_H

#include "array1.h"
#include <algorithm>
#include <cassert>
#include <vector>

template<class T, class ArrayT=std::vector<T> >
struct Array3
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

   int ni, nj, nk;
   ArrayT a;

   // the interface

   Array3(void)
	  : ni(0), nj(0), nk(0)
   {}

   Array3(int ni_, int nj_, int nk_)
	  : ni(ni_), nj(nj_), nk(nk_), a(ni_*nj_*nk_)
   { assert(ni_>=0 && nj_>=0 && nk_>=0); }

   Array3(int ni_, int nj_, int nk_, ArrayT& a_)
	  : ni(ni_), nj(nj_), nk(nk_), a(a_)
   { assert(ni_>=0 && nj_>=0 && nk_>=0); }

   Array3(int ni_, int nj_, int nk_, const T& value)
	  : ni(ni_), nj(nj_), nk(nk_), a(ni_*nj_*nk_, value)
   { assert(ni_>=0 && nj_>=0 && nk_>=0); }

   Array3(int ni_, int nj_, int nk_, const T& value, size_type max_n_)
	  : ni(ni_), nj(nj_), nk(nk_), a(ni_*nj_*nk_, value, max_n_)
   { assert(ni_>=0 && nj_>=0 && nk_>=0); }

   Array3(int ni_, int nj_, int nk_, T* data_)
	  : ni(ni_), nj(nj_), nk(nk_), a(ni_*nj_*nk_, data_)
   { assert(ni_>=0 && nj_>=0 && nk_>=0); }

   Array3(int ni_, int nj_, int nk_, T* data_, size_type max_n_)
	  : ni(ni_), nj(nj_), nk(nk_), a(ni_*nj_*nk_, data_, max_n_)
   { assert(ni_>=0 && nj_>=0 && nk_>=0); }

   template<class OtherArrayT>
   Array3(Array3<T, OtherArrayT>& other)
	   : ni(other.ni), nj(other.nj), nk(other.nk), a(other.a)
   {}


   ~Array3(void)
   {
#ifndef NDEBUG
	  ni=nj=0;
#endif
   }

   // single index reference
   const T& operator()(int index) const
   {
	   assert(index>=0 && index < ni*nj*nk);
	   return a[index];
   }

   T& operator()(int index)
   {
	   assert(index>=0 && index < ni*nj*nk);
	   return a[index];
   }

   const T& operator()(int i, int j, int k) const
   {
	  assert(i>=0 && i<ni && j>=0 && j<nj && k>=0 && k<nk);
	  return a[i+ni*(j+nj*k)];
   }

   T& operator()(int i, int j, int k)
   {
	  assert(i>=0 && i<ni && j>=0 && j<nj && k>=0 && k<nk); 
	  return a[i+ni*(j+nj*k)];
   }

   bool operator==(const Array3<T>& x) const
   { return ni==x.ni && nj==x.nj && nk==x.nk && a==x.a; }

   bool operator!=(const Array3<T>& x) const
   { return ni!=x.ni || nj!=x.nj || nk!=x.nk || a!=x.a; }

   bool operator<(const Array3<T>& x) const
   {
	  if(ni<x.ni) return true; else if(ni>x.ni) return false;
	  if(nj<x.nj) return true; else if(nj>x.nj) return false;
	  if(nk<x.nk) return true; else if(nk>x.nk) return false;
	  return a<x.a;
   }

   bool operator>(const Array3<T>& x) const
   {
	  if(ni>x.ni) return true; else if(ni<x.ni) return false;
	  if(nj>x.nj) return true; else if(nj<x.nj) return false;
	  if(nk>x.nk) return true; else if(nk<x.nk) return false;
	  return a>x.a;
   }

   bool operator<=(const Array3<T>& x) const
   {
	  if(ni<x.ni) return true; else if(ni>x.ni) return false;
	  if(nj<x.nj) return true; else if(nj>x.nj) return false;
	  if(nk<x.nk) return true; else if(nk>x.nk) return false;
	  return a<=x.a;
   }

   bool operator>=(const Array3<T>& x) const
   {
	  if(ni>x.ni) return true; else if(ni<x.ni) return false;
	  if(nj>x.nj) return true; else if(nj<x.nj) return false;
	  if(nk>x.nk) return true; else if(nk<x.nk) return false;
	  return a>=x.a;
   }

   void assign(const T& value)
   { a.assign(value); }

   void assign(int ni_, int nj_, int nk_, const T& value)
   {
	  a.assign(ni_*nj_*nk_, value);
	  ni=ni_;
	  nj=nj_;
	  nk=nk_;
   }
	
   void assign(int ni_, int nj_, int nk_, const T* copydata)
   {
	  a.assign(ni_*nj_*nk_, copydata);
	  ni=ni_;
	  nj=nj_;
	  nk=nk_;
   }
	
   const T& at(int i, int j, int k) const
   {
	  assert(i>=0 && i<ni && j>=0 && j<nj && k>=0 && k<nk);
	  return a[i+ni*(j+nj*k)];
   }

   T& at(int i, int j, int k)
   {
	  assert(i>=0 && i<ni && j>=0 && j<nj && k>=0 && k<nk);
	  return a[i+ni*(j+nj*k)];
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
	  ni=nj=nk=0;
   }

   bool empty(void) const
   { return a.empty(); }

   const_iterator end(void) const
   { return a.end(); }

   iterator end(void)
   { return a.end(); }

   void fill(int ni_, int nj_, int nk_, const T& value)
   {
	  a.fill(ni_*nj_*nk_, value);
	  ni=ni_;
	  nj=nj_;
	  nk=nk_;
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

   void reserve(int reserve_ni, int reserve_nj, int reserve_nk)
   { a.reserve(reserve_ni*reserve_nj*reserve_nk); }

   void resize(int ni_, int nj_, int nk_)
   {
	  assert(ni_>=0 && nj_>=0 && nk_>=0);
	  a.resize(ni_*nj_*nk_);
	  ni=ni_;
	  nj=nj_;
	  nk=nk_;
   }

   void resize(int ni_, int nj_, int nk_, const T& value)
   {
	  assert(ni_>=0 && nj_>=0 && nk_>=0);
	  a.resize(ni_*nj_*nk_, value);
	  ni=ni_;
	  nj=nj_;
	  nk=nk_;
   }

   void set_zero(void)
   { a.set_zero(); }

   size_type size(void) const
   { return a.size(); }

   void swap(Array3<T>& x)
   {
	  std::swap(ni, x.ni);
	  std::swap(nj, x.nj);
	  std::swap(nk, x.nk);
	  a.swap(x.a);
   }

   void trim(void)
   { a.trim(); }

   T trilerp(int i, int j, int k, T fx, T fy, T fz) const {
	   
	   return
		/*  ::trilerp(
		   (*this)(i,j,k), (*this)(i+1,j,k), (*this)(i,j+1,k), (*this)(i+1,j+1,k), 
		   (*this)(i,j,k+1), (*this)(i+1,j,k+1), (*this)(i,j+1,k+1), (*this)(i+1,j+1,k+1), 
		   fx,fy,fz);*/
		 (((*this)(i,j,k) * (1-fx) + (*this)(i+1,j,k) * fx) * (1-fy) + ((*this)(i,j+1,k) * (1-fx) + (*this)(i+1,j+1,k) * fx) * (fy)) * (1-fz) +
		 (((*this)(i,j,k+1) * (1-fx) + (*this)(i+1,j,k+1) * fx) * (1-fy) + ((*this)(i,j+1,k+1) * (1-fx) + (*this)(i+1,j+1,k+1) * fx) * (fy)) * (fz);
   }

   T monotonic_cubic_interp(int i, int j, int k, T fx, T fy, T fz) {

	   if(i<=0||i>=ni-2||j<=0||j>=nj-2||k<=0||k>=nk-2)
		   return 
			 (((*this)(i,j,k) * (1-fx) + (*this)(i+1,j,k) * fx) * (1-fy) + ((*this)(i,j+1,k) * (1-fx) + (*this)(i+1,j+1,k) * fx) * (fy)) * (1-fz) +
			 (((*this)(i,j,k+1) * (1-fx) + (*this)(i+1,j,k+1) * fx) * (1-fy) + ((*this)(i,j+1,k+1) * (1-fx) + (*this)(i+1,j+1,k+1) * fx) * (fy)) * (fz);

	   int ineg1 = i-1; if(ineg1<0) ineg1 = 0;
	   int iplus1 = i+1; if(iplus1 >= ni) iplus1 = ni-1;
	   int iplus2 = i+2; if(iplus2 >= ni) iplus2 = ni-1;
	   int jneg1 = j-1; if(jneg1<0) jneg1 = 0;
	   int jplus1 = j+1; if(jplus1 >= nj) jplus1 = nj-1;
	   int jplus2 = j+2; if(jplus2 >= nj) jplus2 = nj-1;
	   int kneg1 = k-1; if(kneg1<0) kneg1 = 0;
	   int kplus1 = k+1; if(kplus1 >= nk) kplus1 = nk-1;
	   int kplus2 = k+2; if(kplus2 >= nk) kplus2 = nk-1;

	   T z[4];
	   T y[4];

	   y[0] = ::monotonic_cubic_interp((*this)(ineg1,jneg1,kneg1), (*this)(i,jneg1,kneg1), (*this)(iplus1,jneg1,kneg1), (*this)(iplus2,jneg1,kneg1), fx);
	   y[1] = ::monotonic_cubic_interp((*this)(ineg1,j,kneg1), (*this)(i,j,kneg1), (*this)(iplus1,j,kneg1), (*this)(iplus2,j,kneg1), fx);
	   y[2] = ::monotonic_cubic_interp((*this)(ineg1,jplus1,kneg1), (*this)(i,jplus1,kneg1), (*this)(iplus1,jplus1,kneg1), (*this)(iplus2,jplus1,kneg1), fx);
	   y[3] = ::monotonic_cubic_interp((*this)(ineg1,jplus2,kneg1), (*this)(i,jplus2,kneg1), (*this)(iplus1,jplus2,kneg1), (*this)(iplus2,jplus2,kneg1), fx);

	   z[0] = ::monotonic_cubic_interp(y[0],y[1],y[2],y[3], fy);

	   y[0] = ::monotonic_cubic_interp((*this)(ineg1,jneg1,k), (*this)(i,jneg1,k), (*this)(iplus1,jneg1,k), (*this)(iplus2,jneg1,k), fx);
	   y[1] = ::monotonic_cubic_interp((*this)(ineg1,j,k), (*this)(i,j,k), (*this)(iplus1,j,k), (*this)(iplus2,j,k), fx);
	   y[2] = ::monotonic_cubic_interp((*this)(ineg1,jplus1,k), (*this)(i,jplus1,k), (*this)(iplus1,jplus1,k), (*this)(iplus2,jplus1,k), fx);
	   y[3] = ::monotonic_cubic_interp((*this)(ineg1,jplus2,k), (*this)(i,jplus2,k), (*this)(iplus1,jplus2,k), (*this)(iplus2,jplus2,k), fx);

	   z[1] = ::monotonic_cubic_interp(y[0],y[1],y[2],y[3], fy);

	   y[0] = ::monotonic_cubic_interp((*this)(ineg1,jneg1,kplus1), (*this)(i,jneg1,kplus1), (*this)(iplus1,jneg1,kplus1), (*this)(iplus2,jneg1,kplus1), fx);
	   y[1] = ::monotonic_cubic_interp((*this)(ineg1,j,kplus1), (*this)(i,j,kplus1), (*this)(iplus1,j,kplus1), (*this)(iplus2,j,kplus1), fx);
	   y[2] = ::monotonic_cubic_interp((*this)(ineg1,jplus1,kplus1), (*this)(i,jplus1,kplus1), (*this)(iplus1,jplus1,kplus1), (*this)(iplus2,jplus1,kplus1), fx);
	   y[3] = ::monotonic_cubic_interp((*this)(ineg1,jplus2,kplus1), (*this)(i,jplus2,kplus1), (*this)(iplus1,jplus2,kplus1), (*this)(iplus2,jplus2,kplus1), fx);

	   z[2] = ::monotonic_cubic_interp(y[0],y[1],y[2],y[3], fy);

	   y[0] = ::monotonic_cubic_interp((*this)(ineg1,jneg1,kplus2), (*this)(i,jneg1,kplus2), (*this)(iplus1,jneg1,kplus2), (*this)(iplus2,jneg1,kplus2), fx);
	   y[1] = ::monotonic_cubic_interp((*this)(ineg1,j,kplus2), (*this)(i,j,kplus2), (*this)(iplus1,j,kplus2), (*this)(iplus2,j,kplus2), fx);
	   y[2] = ::monotonic_cubic_interp((*this)(ineg1,jplus1,kplus2), (*this)(i,jplus1,kplus2), (*this)(iplus1,jplus1,kplus2), (*this)(iplus2,jplus1,kplus2), fx);
	   y[3] = ::monotonic_cubic_interp((*this)(ineg1,jplus2,kplus2), (*this)(i,jplus2,kplus2), (*this)(iplus1,jplus2,kplus2), (*this)(iplus2,jplus2,kplus2), fx);

	   z[3] = ::monotonic_cubic_interp(y[0],y[1],y[2],y[3], fy);

	   return ::monotonic_cubic_interp(z[0],z[1],z[2],z[3], fz);
   }

   void copy_to(Array3 &other) const {
	   //std::memcpy(a.a, a, ni*nj*sizeof(T));
	   assert(ni == other.ni && nj == other.nj && nk == other.nk);
	  // other.a = a;
	   for(int i = 0; i < ni*nj*nk; i++)
		   other.a[i] = a[i];
   }

   T infnorm() const {
	   T r = 0;
	   for(int i = 0; i < ni*nj*nk; i++) {
		   if(std::fabs(a[i]) > r)
			   r = std::fabs(a[i]);
	   }
	   return r;
   }

};

// some common arrays

typedef Array3<double, Array1<double> > Array3d;
typedef Array3<float, Array1<float> > Array3f;
typedef Array3<long long, Array1<long long> > Array3ll;
typedef Array3<unsigned long long, Array1<unsigned long long> > Array3ull;
typedef Array3<int, Array1<int> > Array3i;
typedef Array3<unsigned int, Array1<unsigned int> > Array3ui;
typedef Array3<short, Array1<short> > Array3s;
typedef Array3<unsigned short, Array1<unsigned short> > Array3us;
typedef Array3<char, Array1<char> > Array3c;
typedef Array3<unsigned char, Array1<unsigned char> > Array3uc;

#endif	// ARRAY3_H
