#include "elastic_wave2D.hpp"
#include "parameters.hpp"

using namespace std;



int main(int argc, char *argv[])
{
  if (argc == 1) // no arguments
  {
    cout << "\nGet help with:\n" << argv[0] << " -h\n" << endl;
    return 0;
  }

  try
  {
    Parameters param;
    param.init(argc, argv);

    ElasticWave2D elwave(param);
    elwave.run();
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



void compute_rayleigh_damping_weights(int n_elements_x,
                                      int n_elements_y,
                                      double left_boundary,
                                      double right_boundary,
                                      double bottom_boundary,
                                      double top_boundary,
                                      double abs_layer_width,
                                      bool left, bool right,
                                      bool bottom, bool top,
                                      double *damping_weights)
{
  const double p = 1.2;

  const double hx = (right_boundary - left_boundary) / n_elements_x;
  const double hy = (top_boundary - bottom_boundary) / n_elements_y;

  for (int ely = 0; ely < n_elements_y; ++ely)
  {
    const double y = bottom_boundary + (ely+0.5)*hy; // center of a cell
    for (int elx = 0; elx < n_elements_x; ++elx)
    {
      const double x = left_boundary + (elx+0.5)*hx; // center of a cell

      double weight = 1e-8; // don't take 0, because of sparse matrix pattern

      if (left && x - abs_layer_width < left_boundary)
        weight += pow((left_boundary-(x-abs_layer_width))/abs_layer_width, p);
      else if (right && x + abs_layer_width > right_boundary)
        weight += pow((x+abs_layer_width-right_boundary)/abs_layer_width, p);

      if (bottom && y - abs_layer_width < bottom_boundary)
        weight += pow((bottom_boundary-(y-abs_layer_width))/abs_layer_width, p);
      else if (top && y + abs_layer_width > top_boundary)
        weight += pow((y+abs_layer_width-top_boundary)/abs_layer_width, p);

      const int el = ely*n_elements_x + elx;
      damping_weights[el] = weight;
    }
  }
}







void get_damp_alphas(double source_frequency, double &alpha1, double &alpha2)
{
  // These numbers have been obtained by Shubin Fu and Kai Gao, while they've
  // been PhD students at Texas A&M.
  const double w1 = 0.7 * source_frequency;
  const double w2 = 20. * source_frequency;
  const double xi1 = 4.5;
  const double xi2 = 0.6;

  alpha1 = 2.*w1*w2*(xi2*w1-xi1*w2)/(w1*w1-w2*w2);
  alpha2 = 2.*(xi1*w1-xi2*w2)/(w1*w1-w2*w2);
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



