#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "Types.h"

class Optimizer
{
public:
  Optimizer(const std::vector<DiffusionEncodingDirection> &encodings);

  void SetDWI(const DiffusionImageType::PixelType &dwi);

  double solve(vnl_vector_fixed<double, 21> &X);

private:
  Optimizer(){};
  unsigned int violated_constraints(vnl_vector_fixed<double, 21> &X);
  double ULLS(vnl_vector_fixed<double, 21> & X);
  double InteriorPointMethod(vnl_vector_fixed<double, 21> & X);
  void Newton(double t, vnl_vector_fixed<double, 21> & X);
  double SolutionNorm(vnl_vector_fixed<double, 21> & X);

  std::vector<DiffusionEncodingDirection> full_encodings;
  std::vector<DiffusionEncodingDirection> nonzero_encodings;
  double bmax;

  static const double t0 = 1e-8;
  static const double tmu = 1.5;
  static const double teps = 0.01;

  vnl_matrix<double> A, C, ATA;
  vnl_vector<double> B, ATB;

  vnl_vector_fixed<double, 21> initialX;
};
