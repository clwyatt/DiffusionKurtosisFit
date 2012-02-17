
#include <vnl/vnl_vector.h>
#include <itkVariableLengthVector.h>

#include "lookup.h"

double tensor_3x3_map(const vnl_vector<double> &x, unsigned int i, unsigned int j)
{
  return x[D[i][j]];
}

double tensor_3x3x3x3_map(const vnl_vector<double> &x, unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
  return x[K[i][j][k][l]];
}

double tensor_3x3_map(const itk::VariableLengthVector<double> &x, unsigned int i, unsigned int j)
{
  return x[D[i][j]];
}

double tensor_3x3x3x3_map(const itk::VariableLengthVector<double> &x, unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
  return x[K[i][j][k][l]];
}

