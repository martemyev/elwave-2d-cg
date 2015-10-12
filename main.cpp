#include "config.hpp"
#include "elastic_wave2D.hpp"
#include "mfem.hpp"
#include "parameters.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;

void show_damp_weights(const Parameters& param);



int main(int argc, char *argv[])
{
  if (argc == 1) // no arguments
  {
    cout << "\nGet help with:\n" << argv[0] << " -h\n" << endl;
    return 0;
  }

#if defined(DEBUG_WAVE)
  cout << "****************************\n";
  cout << "*     DEBUG VERSION        *\n";
  cout << "****************************\n";
#endif

  try
  {
    StopWatch chrono;
    chrono.Start();

    Parameters param;
    param.init(argc, argv);

//    show_SRM_damp_weights(param);

    ElasticWave2D elwave(param);
    elwave.run();

    cout << "\nTOTAL TIME " << chrono.RealTime() << " sec\n" << endl;
  }
  catch (int ierr)
  {
    return ierr;
  }
  catch (...)
  {
    cerr << "\nEXCEPTION\n";
  }

  return 0;
}

/*
void test_polynomials()
{
  int order = 5;

  double *nodes = new double[order+1];

  Poly_1D::UniformPoints(order, nodes);
//  Poly_1D::GaussLobattoPoints(order, nodes);
//  Poly_1D::ChebyshevPoints(order, nodes);
//  Poly_1D::GaussPoints(order, nodes);
  for (int i = 0; i < order+1; ++i)
    cout << nodes[i] << " ";
  cout << endl;

  Poly_1D::Basis *_poly_1D_basis;

  int mode = 0; // 1 means always Lagrange, 0 - Legendre
  if (mode)
    _poly_1D_basis = new Poly_1D::Basis(order, nodes);

  int n_points = 101; //order+1; // 101;
  double *coords = new double[n_points];
//  Poly_1D::UniformPoints(n_points-1, coords);
  Poly_1D::GaussLobattoPoints(n_points-1, coords);
//  Poly_1D::ChebyshevPoints(n_points-1, coords);
//  Poly_1D::GaussPoints(n_points-1, coords);

  Vector *basis = new Vector[n_points];
  cout << "coord=[";
  for (int i = 0; i < n_points; ++i)
  {
    cout << coords[i] << (i == n_points-1 ? "];\n" : ",");
    basis[i].SetSize(order+1);
    if (mode)
      _poly_1D_basis->Eval(coords[i], basis[i]);
    else
      Poly_1D::CalcBasis(order, coords[i], basis[i], Poly_1D::Legendre);
  }

  for (int p = 0; p < order+1; ++p)
  {
    cout << "basis" << p << "=[";
    for (int i = 0; i < n_points; ++i)
      cout << basis[i](p) << (i == n_points-1 ? "];\n" : ",");
  }

//  for (int i = 0; i < n_points; ++i)
//  {
//    cout << "point " << i << ": " << coords[i] << endl;
    for (int p1 = 0; p1 < order; ++p1)
    {
      for (int p2 = 0; p2 < order; ++p2)
        cout << L2(basis, n_points, p1, p2) << " ";
      cout << endl;
    }
//      cout << "L2(basis[" << p << "],basis[" << p+1 << "])=" << L2(basis, n_points, p, p+1) << endl;
//  }

  delete[] basis;
  delete[] coords;
  if (mode)
    delete _poly_1D_basis;
  delete[] nodes;
}
*/


