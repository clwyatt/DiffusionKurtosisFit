#ifndef _Types_h_
#define _Types_h_

#include <itkImage.h>
#include <itkVectorImage.h>

typedef short PixelType;
const unsigned int DIMENSION = 3;
typedef itk::VectorImage<PixelType,DIMENSION> DiffusionImageType;

typedef float InternalPixelType;
typedef itk::Image<InternalPixelType,DIMENSION> B0ImageType;
typedef itk::Image<InternalPixelType,DIMENSION> MDImageType;

typedef itk::VectorImage<InternalPixelType,DIMENSION> TensorImageType;

struct DiffusionEncodingDirection
{
  double bvalue;
  double gx, gy, gz;
  static const double EPS = 0.001;

  bool isZero()
    {
      return ( (fabs(bvalue) < EPS) || ((fabs(gx) < EPS) && (fabs(gy) < EPS) && (fabs(gz) < EPS)) );
    }
};

#endif
