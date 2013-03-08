#include <iostream>

#include "Optimizer.h"

#include "make_A_func.h"
#include "make_C_func.h"

double norm(vnl_vector<double> & x)
{
  double n = 0;
  for(unsigned int i = 0; i < x.size(); ++i) n += x[i]*x[i];
  return sqrt(n);
}

Optimizer::Optimizer(const std::vector<DiffusionEncodingDirection> &encodings)
{
  full_encodings = encodings;
  bmax = 0;
  for(unsigned int i = 0; i < full_encodings.size(); ++i)
    {
      DiffusionEncodingDirection n = full_encodings[i];
      if( !n.isZero() ) nonzero_encodings.push_back(n);
      if(n.bvalue > bmax) bmax = n.bvalue;
    }

  unsigned int N = nonzero_encodings.size();

  // A is N x 21
  A = vnl_matrix<double>(N, 21);
  for(unsigned int i = 0; i < N; ++i)
    {
      DiffusionEncodingDirection n = nonzero_encodings[i];
      make_A(i, n, A);
    }

  // precompute some terms
  ATA = vnl_transpose(A)*A;

  // C is 3N x 21
  C = vnl_matrix<double>(3*N, 21);
  unsigned int r = 0;
  for(unsigned int i = 0; i < N; ++i)
    {
      make_C_block1(r, nonzero_encodings[i], C);
      r += 1;
    }
  for(unsigned int i = 0; i < N; ++i)
    {
      make_C_block2(r, nonzero_encodings[i], C);
      r += 1;
    }
  for(unsigned int i = 0; i < N; ++i)
    {
      make_C_block3(r, nonzero_encodings[i], bmax, C);
      r += 1;
    }

  // B is N x 1
  B = vnl_vector<double>(N);

}

void Optimizer::SetDWI(double * data)
{
  double B0 = 0.;
  unsigned int numZeroEncodings = 0;
  for(unsigned int i = 0; i < full_encodings.size(); ++i)
    {
      DiffusionEncodingDirection n = full_encodings[i];
      if(n.isZero())
	{
	  numZeroEncodings += 1;
	  B0 += data[i];
	}
    }

  B0 = B0/numZeroEncodings;

  unsigned int d = 0;
  for(unsigned int i = 0; i < full_encodings.size(); ++i)
    {
      DiffusionEncodingDirection n = full_encodings[i];
      if(!n.isZero())
	{
	  B[d] = log(data[i]/B0);
	  d += 1;
	}
    }

  // precompute
  ATB = vnl_transpose(A)*B;
}

unsigned int Optimizer::violated_constraints(vnl_vector_fixed<double, 21> &X)
{
  vnl_vector<double> d = C*X;
  unsigned int active_constraints = 0;
  for(unsigned int i = 0; i < d.size(); ++i)
    {
      if(d[i] > 0) active_constraints += 1;
    }
  return active_constraints;
}

double Optimizer::ULLS(vnl_vector_fixed<double, 21> & X)
{
  unsigned int N = nonzero_encodings.size();
  vnl_vector<double> R(N);

  X = vnl_svd<double>(A).solve(B);

  R = A*X-B;

  double norm = 0;
  for(unsigned int i = 0; i < N; ++i)
    {
      norm += R[i]*R[i];
    }
  norm = sqrt(norm);

  return norm;
}

void Optimizer::Newton(double t, vnl_vector_fixed<double, 21> & X)
{
  // std::cout << "Starting Newton with t = " << t;
  // std::cout << " and solution norm = " << SolutionNorm(X) << std::endl;

  unsigned int N = nonzero_encodings.size();

  vnl_vector<double> grad;
  vnl_matrix<double> H;
  vnl_matrix<double> Hinit = vnl_diag_matrix<double>(21,2*t)*ATA;
  double updatenorm;
  const double eps = 1e-6;
  do
    {
      vnl_vector<double> d = C*X;
      grad = vnl_diag_matrix<double>(21,2*t)*(ATA*X-ATB);
      for(unsigned int i = 0; i < 3*N; ++i)
	{
	  grad += (C.get_row(i))/(-d[i]);
	}

      H = Hinit;
      for(unsigned int i = 0; i < 3*N; ++i)
	{
	  vnl_vector<double> temp = C.get_row(i);
	  double constraintsq = d[i]*d[i];
	  for(unsigned int r = 0; r < 21; ++r)
	    for(unsigned int c = 0; c < 21; ++c)
	      {
		H[r][c] += (temp[r]*temp[c])/(constraintsq);
	      }
	}

      vnl_vector<double> update =  vnl_svd<double>(H).solve(grad);
      X = X - update;

      updatenorm = norm(update);

      //std::cout << "Newton update norm: " << updatenorm << std::endl;
    }
  while(updatenorm > eps);

}

double Optimizer::SolutionNorm(vnl_vector_fixed<double, 21> & X)
{
  unsigned int N = nonzero_encodings.size();

  vnl_vector<double> R = A*X-B;
  double norm = 0;
  for(unsigned int i = 0; i < N; ++i)
    {
      norm += R[i]*R[i];
    }

  return sqrt(norm);
}

double Optimizer::InteriorPointMethod(vnl_vector_fixed<double, 21> & X)
{
  unsigned int N = nonzero_encodings.size();

  double t = t0;
  while( (3*N)/t > teps )
    {
      Newton(t, X);
      t = tmu*t;
      if(violated_constraints(X) > 0)
	{
	  std::cout << "WARNING: Interior Point Method Halted Early." << std::endl;
	  break;
	}
    }

  return SolutionNorm(X);
}

double Optimizer::solve(vnl_vector_fixed<double, 21> &X)
{
  vnl_vector_fixed<double, 21> T;
  double norm = ULLS(T);
  if(violated_constraints(T) > 0)
    {
      norm = InteriorPointMethod(X);
    }
  else
    {
      X = T;
    }

  return norm;
}
