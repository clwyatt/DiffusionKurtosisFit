#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "Types.h"

class Optimizer
{
public:
  Optimizer(const std::vector<DiffusionEncodingDirection> &encodings);
  ~Optimizer();

  void SetDWI(double * data);

  double solve(vnl_vector_fixed<double, 21> &X, bool & badfit);

private:
  Optimizer(){};
  void build_ATA();
  void fill_ATB();
  void build_CCT();
  void build_d();
  unsigned int violated_constraints(vnl_vector_fixed<double, 21> &X);
  double ULLS(vnl_vector_fixed<double, 21> & X);
  double InteriorPointMethod(vnl_vector_fixed<double, 21> & X, bool & badfit);
  void Newton(unsigned int timeIndex, vnl_vector_fixed<double, 21> & X);
  double SolutionNorm(vnl_vector_fixed<double, 21> & X);

  std::vector<DiffusionEncodingDirection> full_encodings;
  std::vector<DiffusionEncodingDirection> nonzero_encodings;
  double bmax;

  double t0;
  double tmu;
  double teps;
  std::vector<double> timeSchedule;

  vnl_matrix<double> A, C;
  vnl_vector<double> B;

  vnl_vector<double> grad;
  vnl_matrix<double> H;

  double * ATA;
  double * ATB;
  double * CCT;
  double * d;
};
