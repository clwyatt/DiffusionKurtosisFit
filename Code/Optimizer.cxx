#include <cmath>
#include <cassert>

#include "Optimizer.h"

// lookup.h defines index maps from symmetric tensors
// to optimization free variable x
// constructed using make_lookup.py
#include "lookup.h"

void optimizer_cost_function::f(vnl_vector< double > const &x, vnl_vector< double > &fx)
{
  assert(x.size() == 21);
  unsigned int lenfx = fx.size();

  // TODO unroll inner loops over tensor indices?
  for(unsigned int r = 0; r < lenfx; ++r)
    {
    double b = (*m_encodings)[r].bvalue;
    double g[3];
    g[0] = (*m_encodings)[r].gx;
    g[1] = (*m_encodings)[r].gy;
    g[2] = (*m_encodings)[r].gz;
    double sum1 = 0;
    for(unsigned int i = 0; i < 3; ++i)
      for(unsigned int j = 0; j < 3; ++j)
	{
	sum1 += g[i]*g[j]*x[D[i][j]];
	}

    double sum2 = 0;
    for(unsigned int i = 0; i < 3; ++i)
      for(unsigned int j = 0; j < 3; ++j)
	for(unsigned int k = 0; k < 3; ++k)
	  for(unsigned int l = 0; l < 3; ++l)
	    {
	    sum2 += g[i]*g[j]*g[k]*g[l]*x[K[i][j][k][l]];
	    }

    fx[r] = log((*m_lhs)[r]) - m_logB0 + b*sum1 - b*b*sum2/6;
    }
}

void optimizer_init_x(vnl_vector<double> & x)
{
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
      {
      if(i == j) x[D[i][j]] = 1;
      else x[D[i][j]] = 0;
      }

  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
      for(unsigned int k = 0; k < 3; ++k)
	for(unsigned int l = 0; l < 3; ++l)
	  {
	  if( (i == j) && (i == k) && (i == l) ) x[K[i][j][k][l]] = 1;
	  else x[K[i][j][k][l]] = 0;
	  }
}

