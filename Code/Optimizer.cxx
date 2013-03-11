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

  // create a timeSchedule for interior-point method
  double t = t0;
  while( (3*N)/t > teps )
    {
      timeSchedule.push_back(t);
      t = tmu*t;
    }

  // A is N x 21
  A = vnl_matrix<double>(N, 21);
  for(unsigned int i = 0; i < N; ++i)
    {
      DiffusionEncodingDirection n = nonzero_encodings[i];
      make_A(i, n, A);
    }

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

  grad = vnl_vector<double>(21);
  H = vnl_matrix<double>(21,21);
  build_ATA();
  build_CCT();
  build_d();
}

Optimizer::~Optimizer()
{
  delete [] ATA;
  delete [] ATB;
  delete [] CCT;
  delete [] d;
}

void Optimizer::build_ATA()
{
  unsigned int T = timeSchedule.size();

  ATA = new double[21*21*T];

  for(unsigned int t = 0; t < T; ++t)
    {
      vnl_matrix<double> TEMP = vnl_diag_matrix<double>(21, 2*timeSchedule[t])*(vnl_transpose(A)*A);
      for(unsigned int r = 0; r < 21; ++r)
	for(unsigned int c = 0; c < 21; ++c)
	  {
	    ATA[t*21*21 + r*21 + c] =  TEMP[r][c];
	  }
    }
}

void Optimizer::build_ATB()
{
  unsigned int T = timeSchedule.size();

  ATB = new double[21*T];

  for(unsigned int t = 0; t < T; ++t)
    {
      vnl_vector<double> TEMP = vnl_diag_matrix<double>(21, 2*timeSchedule[t])*(vnl_transpose(A)*B);
      for(unsigned int c = 0; c < 21; ++c)
	  {
	    ATB[t*21 + c] =  TEMP[c];
	  }
    }
}

void Optimizer::build_CCT()
{
  unsigned int N = nonzero_encodings.size();

  CCT = new double[21*21*3*N];

  for(unsigned int i = 0; i < 3*N; ++i)
    {
      vnl_matrix<double> TEMP = outer_product(C.get_row(i), C.get_row(i));
      for(unsigned int r = 0; r < 21; ++r)
	for(unsigned int c = 0; c < 21; ++c)
	  {
	    CCT[i*21*21 + r*21 + c] =  TEMP[r][c];
	  }
    }
}

void Optimizer::build_d()
{
  unsigned int N = nonzero_encodings.size();
  d = new double[3*N];
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
  build_ATB();
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

void Optimizer::Newton(unsigned int timeIndex, vnl_vector_fixed<double, 21> & X)
{
  unsigned int N = nonzero_encodings.size();

  double updatenorm;
  const double eps = 1e-6;
  do
    {
      for(unsigned int r = 0; r < 3*N; ++r)
	{
	  d[r] = 0;
	  for(unsigned int c = 0; c < 21; ++c) d[r] += C[r][c]*X[c];
	}

      for(unsigned int i = 0; i < 21; ++i)
	{
	  grad[i] = 0;
	  for(unsigned int j = 0; j < 21; ++j)
	    {
	      grad[i] += ATA[21*21*timeIndex + 21*i + j]*X[j];
	    }
	  grad[i] -= ATB[21*timeIndex + i];
	}

      for(unsigned int i = 0; i < 3*N; ++i)
	{
	  for(unsigned int j = 0; j < 21; ++j)
	    grad[j] += -C[i][j]/d[i];
	}

      for(unsigned int r = 0; r < 21; ++r)
	for(unsigned int c = 0; c < 21; ++c)
	  {
	    H[r][c] = ATA[21*21*timeIndex + 21*r + c];
	  }

      for(unsigned int i = 0; i < 3*N; ++i)
	{
	  double constraintsq = d[i]*d[i];
	  for(unsigned int r = 0; r < 21; ++r)
	    for(unsigned int c = 0; c < 21; ++c)
	      {
		H[r][c] += CCT[i*21*21 + r*21 + c]/constraintsq;
	      }
	}

      vnl_vector<double> update = vnl_svd<double>(H).solve(grad);
      X = X - update;
      updatenorm = norm(update);
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

double Optimizer::InteriorPointMethod(vnl_vector_fixed<double, 21> & X, bool & badfit)
{
  for(unsigned int timeIndex = 0; timeIndex < timeSchedule.size(); ++timeIndex)
    {
      Newton(timeIndex, X);
      if(violated_constraints(X) > 0)
	{
	  //std::cout << "WARNING: Interior Point Method Halted Early." << std::endl;
	  badfit = true;
	  break;
	}
    }

  return SolutionNorm(X);
}

double Optimizer::solve(vnl_vector_fixed<double, 21> &X, bool & badfit)
{
  badfit = false;

  vnl_vector_fixed<double, 21> T;
  double norm = ULLS(T);
  if(violated_constraints(T) > 0)
    {
      norm = InteriorPointMethod(X, badfit);
    }
  else
    {
      X = T;
    }

  return norm;
}
