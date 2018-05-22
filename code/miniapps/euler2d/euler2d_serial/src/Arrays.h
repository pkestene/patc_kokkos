/**
 * \file Arrays.h
 * \brief Provides CPU/GPU C++ array classes.
 *
 * \author F. Chateau and P. Kestener
 *
 */
#ifndef ARRAYS_H_
#define ARRAYS_H_

#include <stdexcept>
#include <iostream>
#include <cstring>

#include "real_type.h"
#include "common_types.h"


/**
 * \class HostArray Arrays.h 
 * \brief Provides an array object with memory allocated on CPU.
 *
 * HostArray is a storage class for 1d/2d/3d vector data. The number of
 * vector component is specified by nvar (number of variables).
 * In the case of finite volume simulation of Euler equations,
 * nvar should return 4 in 2D and 5 in 3D, as we use the following scalar fields : 
 * rho, E, u, v and w (primitive variables).
 */
template<typename T>
class HostArray
{
public:
  HostArray();
  ~HostArray();

  /** memory allocation for 1D data */
  void allocate(int length, int numVar);
  /** memory allocation for 2D data */
  void allocate(uint3 dim);
  /** memory allocation for 3D data */
  void allocate(uint4 dim);
  /** memory free */
  void free();

  /** copy from another array (make a call to allocate and then copy data) */
  void copyHard(HostArray<T>& src);

  /** copy to another existing array of the same size */
  void copyTo(HostArray<T>& src);

  uint dimx() const	{ return _dim.x; }
  uint dimy() const	{ return _dim.y; }
  uint dimz() const	{ return _dim.z; }
  uint nvar() const	{ return _dim.w; }

  uint pitch() const	{ return _dim.x; }
  uint section() const	{ return pitch() * _dim.y * _dim.z; }
  uint size() const	{ return pitch() * _dim.y * _dim.z * _dim.w; }

  uint dimXBytes() const	{ return dimx() * sizeof(T); }
  uint pitchBytes() const 	{ return pitch() * sizeof(T); }
  uint sectionBytes() const	{ return section() * sizeof(T); }
  uint sizeBytes() const	{ return size()  * sizeof(T); }

  T* data()		{ return _data; }
  const T* data() const	{ return _data; }

  /** access 1d data (only valied if _dim.y and _dim.z are 1)*/
  T& operator() (int i, int ivar) { 
    return _data[i+_dim.x*ivar]; }
  T  operator() (int i, int ivar) const { 
    return _data[i+_dim.x*ivar]; }

  /** access 2d data (only valid if _dim.z=1) */
  T& operator() (int i, int j, int ivar) { 
    return _data[i+_dim.x*(j+_dim.y*ivar)]; }
  T  operator() (int i, int j, int ivar) const { 
    return _data[i+_dim.x*(j+_dim.y*ivar)]; }

  /** access 3d data */
  T& operator() (int i, int j, int k, int ivar) { 
    return _data[i+_dim.x*(j+_dim.y*(k+_dim.z*ivar))]; }
  T  operator() (int i, int j, int k, int ivar) const { 
    return _data[i+_dim.x*(j+_dim.y*(k+_dim.z*ivar))]; }

  /** access data directly */
  T& operator() (int i) { return _data[i]; }
  T  operator() (int i) const { return _data[i]; }

  /** other operators */
  HostArray<T> &operator+=(const HostArray<T>& operand);
  HostArray<T> &operator-=(const HostArray<T>& operand);
  HostArray<T> &operator*=(const HostArray<T>& operand);
  HostArray<T> &operator/=(const HostArray<T>& operand);

  /** other methods */
  void reset() {memset((void*) _data, 0, sizeBytes());  };

  /** init with constant value */
  void init(T value) {
    for (int idx=0; idx<this->size(); idx++)
      _data[idx] = value;
  };

private:
  T*	_data;
  uint4	_dim;

public:
  /** is memory allocated ? */
  bool isAllocated;

  /** total allocated memory in bytes */
  static unsigned long int totalAllocMemoryInKB;

}; // class HostArray

template<typename T> unsigned long int HostArray<T>::totalAllocMemoryInKB = 0;

template<typename T> bool arraysHaveSameShape(const HostArray<T> &array1,
					      const HostArray<T> &array2) {
  if (array1.dimx() != array2.dimx() ||
      array1.dimy() != array2.dimy() ||
      array1.dimz() != array2.dimz() ||
      array1.nvar() != array2.nvar() )
    return false;
  return true;
} // arraysHaveSameShape

  // =======================================================
  // =======================================================
  /**
   * A simple routine to print a HostArray to an ostream
   */
#define PRECISION 7
#define WIDTH     10
template<typename T>
std::ostream& operator<<(std::ostream& os, const HostArray<T>& U) 
{
    
  // print a HostArray
  os << "HostArray values:" << std::endl;
  {
    for (uint nVar = 0; nVar < U.nvar(); ++nVar) {
      os << "nVar = " << nVar << std::endl;
      for (uint k = 0; k < U.dimz(); ++k) {
	os << "k = " << k << std::endl;
	for (uint j = 0; j < U.dimy(); ++j) {
	  for (uint i = 0; i < U.dimx(); ++i) {
	    os.precision(PRECISION); os.width(WIDTH);
	    os << static_cast<T>(U(i,j,k,nVar)) << " ";
	  }
	  os << std::endl;
	}
	os << std::endl;
      }
      os << std::endl;
    }
    return os;
  }
} // operator<<




////////////////////////////////////////////////////////////////////////////////
// HostArray class methods body
////////////////////////////////////////////////////////////////////////////////

// =======================================================
// =======================================================
template<typename T>
HostArray<T>::HostArray()
  : _data(0), _dim(make_uint4(0, 0, 0, 0)), isAllocated(false)
{
}

// =======================================================
// =======================================================
template<typename T>
HostArray<T>::~HostArray()
{
  free();
}

// =======================================================
// =======================================================
/* 1d data allocation */
template<typename T>
void HostArray<T>::allocate(int length, 
			    int numVar)
{
    
  free();
  _dim.x = length;
  _dim.y = 1;
  _dim.z = 1;
  _dim.w = numVar;
  _data = new T[length * numVar];
    
  isAllocated = true;

  totalAllocMemoryInKB += (length * numVar * sizeof(T) / 1024);   

} // HostArray<T>::allocate for 1D data
  
  // =======================================================
  // =======================================================
  /* 2d data allocation */
template<typename T>
void HostArray<T>::allocate(uint3 dim)
{

  free();
  _dim.x = dim.x;
  _dim.y = dim.y;
  _dim.z = 1;
  _dim.w = dim.z;
  _data = new T[dim.x * dim.y * dim.z];

  isAllocated = true;

  totalAllocMemoryInKB += (dim.x * dim.y * dim.z * sizeof(T) / 1024);

} // void HostArray<T>::allocate for 2D data

  // =======================================================
  // =======================================================
/* 3d data allocation */
template<typename T>
void HostArray<T>::allocate(uint4 dim)
{

  free();
  _dim = dim;
  _data = new T[dim.x * dim.y * dim.z * dim.w];

  isAllocated = true;

  totalAllocMemoryInKB += (dim.x * dim.y * dim.z * dim.w * sizeof(T) / 1024);

} // void HostArray<T>::allocate for 3D data

  // =======================================================
  // =======================================================
template<typename T>
void HostArray<T>::free()
{

  delete[] _data;

  isAllocated = false;

} // HostArray<T>::free

  // =======================================================
  // =======================================================
template<typename T>
void HostArray<T>::copyHard(HostArray<T>& src)
{
  uint4	src_dim;

  src_dim.x = src.dimx();
  src_dim.y = src.dimy();
  src_dim.z = src.dimz();
  src_dim.w = src.nvar();
  
  // memory allocation (eventually free before allocate, if previously allocated)
  if (src_dim.y == 1 and src_dim.z == 1) { // ONE_D
    this->allocate(src_dim.x, src_dim.w);
  } else if (src_dim.y != 1 and src_dim.z == 1) { // TWO_D
    this->allocate(make_uint3(src_dim.x, src_dim.y, src_dim.w));
  } else { // THREE_D
    this->allocate(src_dim);
  }
  
  // copy data
  //T* src_data = src.data();
  for (unsigned int i=0; i<src.size(); i++) {
    this->_data[i] = src(i);
  }

} // HostArray<T>::copyHard

  // =======================================================
  // =======================================================
template<typename T>
void HostArray<T>::copyTo(HostArray<T>& dest)
{

  if (_dim.x != dest.dimx() or
      _dim.y != dest.dimy() or
      _dim.z != dest.dimz() or
      _dim.w != dest.nvar() ) {
    std::cerr << "HostArray dimensions do not match ! abort...\n";
    return;
  }
  
  // copy data
  //T* dest_data = dest.data();
  for (unsigned int i=0; i<dest.size(); i++) {
    dest(i) = this->_data[i];
  }

} // HostArray<T>::copyTo

  // =======================================================
  // =======================================================
template<typename T>
HostArray<T> &HostArray<T>::operator+=(const HostArray<T> &operand)
{

  // check arrays have same shape
  if (!arraysHaveSameShape(*this, operand)) {
    std::cerr << "HostArray<T>::operator+= : arrays do not have same shape\n";
    return *this;
  }
   
  // apply operator
  for (unsigned int i=0; i<size(); i++) {
    this->_data[i] += operand(i);
  }

  return *this;

} // HostArray<T>::operator+=

  // =======================================================
  // =======================================================
template<typename T>
HostArray<T> &HostArray<T>::operator-=(const HostArray<T> &operand)
{

  // check arrays have same shape
  if (!arraysHaveSameShape(*this, operand)) {
    std::cerr << "HostArray<T>::operator-= : arrays do not have same shape\n";
    return *this;
  }
   
  // apply operator
  for (unsigned int i=0; i<size(); i++) {
    this->_data[i] -= operand(i);
  }

  return *this;

} // HostArray<T>::operator-=

  // =======================================================
  // =======================================================
template<typename T>
HostArray<T> &HostArray<T>::operator*=(const HostArray<T> &operand)
{

  // check arrays have same shape
  if (!arraysHaveSameShape(*this, operand)) {
    std::cerr << "HostArray<T>::operator*= : arrays do not have same shape\n";
    return *this;
  }
   
  // apply operator
  for (unsigned int i=0; i<size(); i++) {
    this->_data[i] *= operand(i);
  }

  return *this;

} // HostArray<T>::operator*=

  // =======================================================
  // =======================================================
template<typename T>
HostArray<T> &HostArray<T>::operator/=(const HostArray<T> &operand)
{

  // check arrays have same shape
  if (!arraysHaveSameShape(*this, operand)) {
    std::cerr << "HostArray<T>::operator/= : arrays do not have same shape\n";
    return *this;
  }
   
  // apply operator
  for (unsigned int i=0; i<size(); i++) {
    this->_data[i] /= operand(i);
  }

  return *this;

} // HostArray<T>::operator/=

#endif /*ARRAYS_H_*/
