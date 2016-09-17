/*=========================================================================

  Module: vtkImageSimilarityMetricInternals.h

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// This is an internal header for use in image similarity metrics.
// It helps the multithreading of the metrics by defining templates
// for thread-local storage, much like vtkSMPTools but also compatible
// with vtkMultiThreader so that either can be used.

#ifndef vtkImageSimilarityMetricInternals_h
#define vtkImageSimilarityMetricInternals_h

#include <vtkThreadedImageAlgorithm.h>

// Do the VTK version check to see if vtkSMPTools will be used
#if VTK_MAJOR_VERSION > 7 || (VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION >= 0)
#define USE_SMP_THREADED_IMAGE_ALGORITHM
#include <vtkSMPTools.h>
#endif

#ifdef USE_SMP_THREADED_IMAGE_ALGORITHM
//----------------------------------------------------------------------------
// Given a per-thread data structure T, this template provides an iterator
// that can be used either with vtkMultiThreader by constructing it with T*,
// or with vtkSMPTools by constructing it with vtkSMPThreadLocal<T>::iterator.
// It allows common code to be written for doing Reduce() operations with
// either vtkSMPTools or with vtkMultiThreader.
template<class T, class TLS>
class vtkImageSimilarityMetricTLSIterator
{
public:
  vtkImageSimilarityMetricTLSIterator(
    const typename TLS::iterator& iter) : Pter(0), Iter(iter) {}

  vtkImageSimilarityMetricTLSIterator(T *iter) : Pter(iter) {}

  vtkImageSimilarityMetricTLSIterator& operator++()
    {
    if (this->Pter) { ++this->Pter; } else { ++this->Iter; }
    return *this;
    }

  vtkImageSimilarityMetricTLSIterator operator++(int)
    {
    vtkImageSimilarityMetricTLSIterator copy = *this;
    if (this->Pter) { ++this->Pter; } else { ++this->Iter; }
    return copy;
    }

  bool operator==(const vtkImageSimilarityMetricTLSIterator& other)
    {
    return (this->Pter ? (this->Pter == other.Pter)
                       : (this->Iter == other.Iter));
    }

  bool operator!=(const vtkImageSimilarityMetricTLSIterator& other)
    {
    return (this->Pter ? (this->Pter != other.Pter)
                       : (this->Iter != other.Iter));
    }

  T& operator*()
    {
    return (this->Pter ? *this->Pter : *this->Iter);
    }

  T* operator->()
    {
    return (this->Pter ? this->Pter : &*this->Iter);
    }

private:
  T *Pter;
  typename TLS::iterator Iter;
};

//----------------------------------------------------------------------------
template<typename T>
class vtkImageSimilarityMetricTLS
{
public:
  typedef vtkImageSimilarityMetricTLSIterator<
    T, vtkSMPThreadLocal<T> > iterator;

  vtkImageSimilarityMetricTLS()
    : MT(0), NumberOfThreads(0) {}

  ~vtkImageSimilarityMetricTLS()
    {
    delete [] this->MT;
    }

  void Initialize(vtkImageSimilarityMetric *a)
    {
    if (!a->GetEnableSMP())
      {
      size_t n = a->GetNumberOfThreads();
      this->MT = new T[n];
      this->NumberOfThreads = n;
      }
    }

  T& Local(size_t threadId)
    {
    if (this->MT == 0)
      {
      return this->SMP.Local();
      }
    else
      {
      return this->MT[threadId];
      }
    }

  iterator begin()
    {
    if (this->MT == 0)
      {
      return this->SMP.begin();
      }
    else
      {
      return this->MT;
      }
    }

  iterator end()
    {
    if (this->MT == 0)
      {
      return this->SMP.end();
      }
    else
      {
      return this->MT + this->NumberOfThreads;
      }
    }

private:
  vtkImageSimilarityMetricTLS(const vtkImageSimilarityMetricTLS&);
  void operator=(const vtkImageSimilarityMetricTLS&);

  vtkSMPThreadLocal<T> SMP;
  T *MT;
  size_t NumberOfThreads;
};

#else

//----------------------------------------------------------------------------
// A simple thread-local storage class template for versions of VTK prior
// to the introduction of vtkSMPTools.  It simply provides an array of the
// specified object type, one object per thread.
template<typename T>
class vtkImageSimilarityMetricTLS
{
public:
  typedef T* iterator;

  vtkImageSimilarityMetricTLS() : MT(0), NumberOfThreads(0) {}

  ~vtkImageSimilarityMetricTLS()
    {
    delete [] this->MT;
    }

  void Initialize(vtkImageSimilarityMetric *a)
    {
    size_t n = a->GetNumberOfThreads();
    this->MT = new T[n];
    this->NumberOfThreads = n;
    }

  T& Local(size_t threadId)
    {
    return this->MT[threadId];
    }

  iterator begin()
    {
    return this->MT;
    }

  iterator end()
    {
    return this->MT + this->NumberOfThreads;
    }

private:
  vtkImageSimilarityMetricTLS(const vtkImageSimilarityMetricTLS&);
  void operator=(const vtkImageSimilarityMetricTLS&);

  T *MT;
  size_t NumberOfThreads;
};

#endif

#endif /* vtkImageSimilarityMetricInternals_h */
