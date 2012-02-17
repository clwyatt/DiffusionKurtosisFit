#include <cmath>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_least_squares_function.h>
#include <vnl/algo/vnl_levenberg_marquardt.h>

#include "Types.h"

class optimizer_cost_function: public vnl_least_squares_function
{
public:

  optimizer_cost_function(vnl_vector<double> *lhs, 
			  double B0, 
			  std::vector<DiffusionEncodingDirection> *encodings)
    : vnl_least_squares_function(21, lhs->size(), vnl_least_squares_function::no_gradient)
    {
      m_lhs = lhs;
      m_encodings = encodings;
      m_logB0 = log(B0);
    }

  void f(vnl_vector< double > const &x, vnl_vector< double > &fx);

private:

  vnl_vector<double> *m_lhs;
  std::vector<DiffusionEncodingDirection> *m_encodings;
  double m_logB0;
};

void optimizer_init_x(vnl_vector<double> & x);
