#include "GLL_quadrature.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace std;
using namespace mfem;

//void segment_GLL_quadrature(int p, double *x, double *w)
//{
//  const int n = p+1; // number of points
//  double *xold = new double[n];

//  // Use the Chebyshev-Gauss-Lobatto nodes as the first guess
//  for (int i = 0; i < n; ++i)
//  {
//    x[i] = cos(M_PI * 1.0*i / (1.0*(n-1)));
//    xold[i] = 2.0;
//  }

//  const double tol = 2e-16;
//  const int maxiter = 10000;

//  double *P = new double[n*n];

//  int iter = 0;
//  while (maxabsdiff(x, xold, n) > tol)
//  {
//    for (int i = 0; i < n; ++i)
//    {
//      xold[i] = x[i];
//      P[0*n+i] = 1.0;
//      P[1*n+i] = x[i];
//    }

//    for (int k = 2; k < n; ++k)
//      for (int i = 0; i < n; ++i)
//        P[k*n+i] = ((2.0*k-1.0)*x[i]*P[(k-1)*n+i] - (k-1.0)*P[(k-2)*n+i]) / k;

//    for (int i = 0; i < n; ++i)
//      x[i] = xold[i] - (x[i]*P[(n-1)*n+i] - P[(n-2)*n+i]) / (n*P[(n-1)*n+i]);

//    if (++iter == maxiter)
//      throw std::runtime_error("segment_GLL_quadrature failed");
//  }

//  delete[] xold;

//  for (int i = 0; i < n; ++i)
//    w[i] = 2.0 / (n*(n-1)*P[(n-1)*n+i]*P[(n-1)*n+i]);

//  delete[] P;
//}

/**
 * ||a-b||_{L^Infinity}
 */
double LInfDiff(const Vector& a, const Vector& b)
{
  Vector diff = a;
  diff -= b;
  return diff.Normlinf();
}

/**
 * Computation of Gauss-Lobatto-Legendre quadrature points at a reference
 * segment [-1, 1] and corresponding weights. Adopted from a Matlab code written
 * by Greg von Winckel.
 * @param p - order
 * @param x - GLL points
 * @param w - GLL weights
 */
void segment_GLL_quadrature(int p, Vector& x, Vector& w, double tol = 2e-16,
                            int maxiter = 10000)
{
  const int n = p+1; // number of points

  x.SetSize(n);
  w.SetSize(n);
  Vector xold(n); // auxiliary vector

  // Use the Chebyshev-Gauss-Lobatto nodes as the first guess
  for (int i = 0; i < n; ++i)
  {
    x[i] = cos(M_PI * 1.0*i / (1.0*(n-1)));
    xold[i] = 2.0;
  }

  double *P = new double[n*n];

  int iter = 0;
  while (LInfDiff(x, xold) > tol)
  {
    for (int i = 0; i < n; ++i)
    {
      xold[i] = x[i];
      P[0*n+i] = 1.0;
      P[1*n+i] = x[i];
    }

    for (int k = 2; k < n; ++k)
      for (int i = 0; i < n; ++i)
        P[k*n+i] = ((2.0*k-1.0)*x[i]*P[(k-1)*n+i] - (k-1.0)*P[(k-2)*n+i]) / k;

    for (int i = 0; i < n; ++i)
      x[i] = xold[i] - (x[i]*P[(n-1)*n+i] - P[(n-2)*n+i]) / (n*P[(n-1)*n+i]);

    MFEM_VERIFY(++iter != maxiter, "Computation of GLL quadrature points "
                "failed");
  }

  for (int i = 0; i < n; ++i)
    w[i] = 2.0 / (n*(n-1)*P[(n-1)*n+i]*P[(n-1)*n+i]);

  delete[] P;
}

void create_segment_GLL_rule(int p, mfem::IntegrationRule& segment_GLL)
{
  segment_GLL.SetSize(p+1);

  { // computation of GLL points in MFEM (for comparison)
    double *GLL_points_mfem= new double[p+1];
    Poly_1D::GaussLobattoPoints(p, GLL_points_mfem);
    cout << "GLL points 2 on [0, 1]:\n";
    for (int i = 0; i < p+1; ++i)
      cout << setw(10) << GLL_points_mfem[i] << endl;
    delete[] GLL_points_mfem;
  }

  Vector GLL_points, GLL_weights;
  segment_GLL_quadrature(p, GLL_points, GLL_weights);

  cout << "GLL points & weights on [-1, 1]:\n";
  for (int i = 0; i < p+1; ++i)
    cout << setw(10) << GLL_points[i] << " " << GLL_weights[i] << endl;

  cout << "GLL points & weights on [0, 1]:\n";
  for (int i = 0; i < p+1; ++i)
  {
    // shift from [-1, 1] to [0, 1]
    segment_GLL.IntPoint(i).x      = 0.5*GLL_points[i] + 0.5;
    segment_GLL.IntPoint(i).weight = 0.5*GLL_weights[i];

    cout << setw(10) << segment_GLL.IntPoint(i).x
         << " " << segment_GLL.IntPoint(i).weight << endl;
  }
}

