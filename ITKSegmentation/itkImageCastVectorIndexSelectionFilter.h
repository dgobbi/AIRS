/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageCastVectorIndexSelectionFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/10/04 20:00:42 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageCastVectorIndexSelectionFilter_h
#define __itkImageCastVectorIndexSelectionFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{
  
namespace Functor {   
  
template< class TInput, class TOutput>
class CastVectorIndexSelection
{
public:
  CastVectorIndexSelection() {m_Index = 0;}
  ~CastVectorIndexSelection() {}

  unsigned int GetIndex() const { return m_Index; }
  void SetIndex(unsigned int i) { m_Index = i; }

  inline TOutput operator()( const TInput & A )
  {
    TOutput B;
    B[m_Index] = A;
    return B;
  }
  
  inline TOutput SetElement( const TInput & A, const TOutput & B )
  {
    return static_cast<TOutput>( B[m_Index] = A );
  }
      
private:
  unsigned int m_Index;   
}; 
}



 /** \class ImageCastVectorIndexSelectionFilter
 *
 * \brief Extracts the selected index of the vector that is the input
 * pixel type
 *
 * This filter is templated over the input image type and 
 * output image type.
 * 
 * The filter expect the input image pixel type to be a vector and 
 * the output image pixel type to be a scalar. The only requirement on
 * the type used for representing the vector is that it must provide an
 * operator[].
 *
 * \ingroup IntensityImageFilters  Multithreaded
 */


template <class TInputImage, class TOutputImage>
class ITK_EXPORT ImageCastVectorIndexSelectionFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Functor::CastVectorIndexSelection< typename TInputImage::PixelType, 
                                                           typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef ImageCastVectorIndexSelectionFilter Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Functor::CastVectorIndexSelection< typename TInputImage::PixelType, 
                                                                     typename TOutputImage::PixelType> > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
    
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Get/Set methods for the index */
  void SetIndex(unsigned int i)
    {
      if (i != this->GetFunctor().GetIndex())
        {
        this->GetFunctor().SetIndex(i);
        this->Modified();
        }
    }
  unsigned int GetIndex(void) const { return this->GetFunctor().GetIndex(); }
  
  typedef typename Superclass::OutputImageRegionType  OutputImageRegionType;
  typedef typename Superclass::InputImagePointer      InputImagePointer;
  typedef typename Superclass::OutputImagePointer     OutputImagePointer;
  typedef typename Superclass::InputImageRegionType   InputImageRegionType;

    
  virtual void ThreadedGenerateData( 
                      const OutputImageRegionType &outputRegionForThread,
                      int threadId);

protected:
  ImageCastVectorIndexSelectionFilter() {}
  virtual ~ImageCastVectorIndexSelectionFilter() {}
    
private:
  ImageCastVectorIndexSelectionFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};


/**
 * ThreadedGenerateData Performs the pixel-wise addition
 */
template <class TInputImage, class TOutputImage >
void
ImageCastVectorIndexSelectionFilter<TInputImage,TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread,
                        int threadId)
{
  InputImagePointer  inputPtr = this->GetInput();
  OutputImagePointer outputPtr = this->GetOutput(0);
  
  // Define the portion of the input to walk for this thread, using
  // the CallCopyOutputRegionToInputRegion method allows for the input
  // and output images to be different dimensions
  InputImageRegionType inputRegionForThread;
  this->CallCopyOutputRegionToInputRegion(inputRegionForThread, outputRegionForThread);

  // Define the iterators
  ImageRegionConstIterator<TInputImage>  inputIt(inputPtr, inputRegionForThread);
  ImageRegionIterator<TOutputImage> outputIt(outputPtr, outputRegionForThread);

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  inputIt.GoToBegin();
  outputIt.GoToBegin();

  typename TOutputImage::PixelType B;
  const unsigned int index = this->GetFunctor().GetIndex();   
  while( !inputIt.IsAtEnd() ) 
    {
    B =  outputIt.Get();
    B[ index ] = inputIt.Get();
    outputIt.Set( B );
    ++inputIt;
    ++outputIt;
    progress.CompletedPixel();  // potential exception thrown here
    }
}

} // end namespace itk


#endif
