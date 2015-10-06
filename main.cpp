#include "elastic_wave2D.hpp"
#include "GLL_quadrature.hpp"
#include "mfem.hpp"
#include "parameters.hpp"
#include "utilities.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;
using namespace mfem;

const double SOURCE_X = 500.0;
const double SOURCE_Y = 500.0;
const double RHO      = 1000.0;
const double VP       = 3000.0;
const double VS       = 1700.0;

const double FREQUENCY     = 20;
const double RICKER_SCALE  = 1e+10;
const double GAUSS_SUPPORT = 1;
const double POINT_FORCE_AMP = 1;
const double SOURCE_ANGLE  = 270;

const double M_XX = 0; //1;
const double M_XY = 0;
const double M_YY = 0; //1;

const double RECEIVER_Y = 700;

const int STEP_SNAP = 500;

void save_vts(int time_step, const std::string &solname,
              double sx, double sy, int nx, int ny,
              const Vector& sol_x, const Vector& sol_y);

void compute_rayleigh_damping_weights(int n_elements_x,
                                      int n_elements_y,
                                      double left_boundary,
                                      double right_boundary,
                                      double bottom_boundary,
                                      double top_boundary,
                                      double abs_layer_width,
                                      bool left, bool right,
                                      bool bottom, bool top,
                                      double *damping_weights);



void get_damp_alphas(double source_frequency, double &alpha1, double &alpha2);

void solve_without_damp_layer(int nx, int ny, double *rho_array,
                              double *vp_array, double *vs_array, int order,
                              int dim, Mesh *mesh, double sx, double sy,
                              double tend, double tstep, bool visualization);

double L2(Vector *basis, int n_points, int p1, int p2)
{
  double sum = 0;
  for (int i = 0; i < n_points; ++i)
    sum += basis[i](p1)*basis[i](p2);
  return sqrt(fabs(sum));
}



int main(int argc, char *argv[])
{
  try
  {
    Parameters param;
    param.init(argc, argv);

    ElasticWave2D elwave(param);
    elwave.run_SEM();
  }
  catch (int ierr)
  {
    return ierr;
  }
  catch (...)
  {
    cerr << "\nEXCEPTION\n";
  }



//  const int nx = 100;
//  const double sx = 1000;
//  const int order = 4;
//  Mesh *mesh = new Mesh(nx, sx);
//  const int dim = mesh->Dimension();
//  cout << "dim = " << dim << endl;

//  IntegrationRule segment_GLL;
//  segment_GLL.SetSize(order+1);

//  double *GLL_points  = new double[order+1];
//  double *GLL_weights = new double[order+1];

//  double *GLL_points_2= new double[order+1];
//  Poly_1D::GaussLobattoPoints(order, GLL_points_2);
//  cout << "GLL points 2 on [0, 1]:\n";
//  for (int i = 0; i < order+1; ++i)
//    cout << setw(10) << GLL_points_2[i] << endl;
//  delete[] GLL_points_2;

//  segment_GLL_quadrature(order, GLL_points, GLL_weights);

//  cout << "GLL points & weights on [-1, 1]:\n";
//  for (int i = 0; i < order+1; ++i)
//    cout << setw(10) << GLL_points[i] << " " << GLL_weights[i] << endl;

//  cout << "GLL points & weights on [0, 1]:\n";
//  for (int i = 0; i < order+1; ++i)
//  {
//    // shift from [-1, 1] to [0, 1]
//    segment_GLL.IntPoint(i).x      = 0.5*GLL_points[i] + 0.5;
//    segment_GLL.IntPoint(i).weight = 0.5*GLL_weights[i];

//    cout << setw(10) << segment_GLL.IntPoint(i).x
//         << " " << segment_GLL.IntPoint(i).weight << endl;
//  }

//  delete[] GLL_weights;
//  delete[] GLL_points;

//  FiniteElementCollection *fec;
//  fec = new H1_FECollection(order, dim);
//  FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec, dim, Ordering::byVDIM);
//  cout << "Number of unknowns: " << fespace->GetVSize() << endl;

//  BilinearForm *m = new BilinearForm(fespace);
//  ConstantCoefficient rho(RHO);
//  MassIntegrator mass_int(rho);
//  mass_int.SetIntRule(&segment_GLL);
//  m->AddDomainIntegrator(&mass_int);
//  m->Assemble();
//  m->Finalize();
//  SparseMatrix& M = m->SpMat();

//  ofstream mout("mass_mat.dat");
//  m->PrintMatlab(mout);

//  cout << "M.nnz = " << M.NumNonZeroElems() << endl;
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
  // These numbers have been obtained by Kai Gao, while he's been a PhD student
  // at Texas A&M.
  const double w1 = 0.7 * source_frequency;
  const double w2 = 20. * source_frequency;
  const double xi1 = 4.5;
  const double xi2 = 0.6;

  alpha1 = 2.*w1*w2*(xi2*w1-xi1*w2)/(w1*w1-w2*w2);
  alpha2 = 2.*(xi1*w1-xi2*w2)/(w1*w1-w2*w2);
}




//void solve_without_damp_layer(int nx, int ny, double *rho_array,
//                              double *vp_array, double *vs_array, int order,
//                              int dim, Mesh *mesh, double sx, double sy,
//                              double tend, double tstep, bool visualization)
//{
//  const int n_elements = nx*ny;
//  double *lambda_array   = new double[n_elements];
//  double *mu_array       = new double[n_elements];
//  for (int i = 0; i < n_elements; ++i)
//  {
//    const double rho = rho_array[i];
//    const double vp  = vp_array[i];
//    const double vs  = vs_array[i];

//    lambda_array[i]  = rho*(vp*vp - 2.*vs*vs);
//    mu_array[i]      = rho*vs*vs;
//  }
//  CWConstCoefficient rho_cw(rho_array);
//  CWConstCoefficient lambda_cw(lambda_array);
//  CWConstCoefficient mu_cw(mu_array);

//  FiniteElementCollection *fec;
//  fec = new H1_FECollection(order, dim);
//  FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec, dim, Ordering::byVDIM);
//  cout << "Number of unknowns: " << fespace->GetVSize() << endl;


//  IntegrationRule segment_GLL;
//  create_segment_GLL_rule(order, segment_GLL);

//  IntegrationRule quad_GLL(segment_GLL, segment_GLL);

//  BilinearForm *s = new BilinearForm(fespace);
//  ElasticityIntegrator elast_int(lambda_cw, mu_cw);
//  elast_int.SetIntRule(&quad_GLL);
//  s->AddDomainIntegrator(&elast_int);
//  s->Assemble();
//  s->Finalize();
//  const SparseMatrix& S = s->SpMat();

//  BilinearForm *m = new BilinearForm(fespace);
//  VectorMassIntegrator mass_int(rho_cw);
//  mass_int.SetIntRule(&quad_GLL);
//  m->AddDomainIntegrator(&mass_int);
//  m->Assemble();
//  m->Finalize();
//  SparseMatrix& M = m->SpMat();

////  ofstream mout("mass_mat.dat");
////  m->PrintMatlab(mout);
////  cout << "M.nnz = " << M.NumNonZeroElems() << endl;

//  Vector diagM;
//  M.GetDiag(diagM); // mass matrix is diagonal

//  Vector invDiagM(diagM.Size());
//  for (int i = 0; i < diagM.Size(); ++i)
//    invDiagM[i] = 1.0 / diagM[i];

//  // Define a simple symmetric Gauss-Seidel preconditioner.
//  GSSmoother prec(M);

//  GridFunction x_0(fespace);
//  GridFunction x_1(fespace);
//  GridFunction x_2(fespace);
//  x_0 = 0.0;
//  x_1 = 0.0;
//  x_2 = 0.0;


//  // Set up the linear form b(.) which corresponds to the RHS.
//  LinearForm *b = new LinearForm(fespace);
//  VectorFunctionCoefficient f(dim, Gauss);
//  VectorDomainLFIntegrator point_source_int(f);
//  point_source_int.SetIntRule(&quad_GLL);
//  b->AddDomainIntegrator(&point_source_int);
//  VectorFunctionCoefficient div_moment_tensor(dim, DivMomentTensor);
//  VectorDomainLFIntegrator dmt_int(div_moment_tensor);
//  dmt_int.SetIntRule(&quad_GLL);
//  b->AddDomainIntegrator(&dmt_int);
//  b->Assemble();

//  cout << "||b||_L2 = " << b->Norml2() << endl;

//  const double hy = sy / ny;
//  const int rec_y_index = RECEIVER_Y / hy;

//  ofstream solout_x("seis0.bin", std::ios::binary);
//  ofstream solout_y("seis1.bin", std::ios::binary);

//  const int n_time_steps = ceil(tend/tstep);
//  const int tenth = 0.1 * n_time_steps;

//  cout << "N time steps = " << n_time_steps
//       << "\nTime loop..." << endl;

//  const int N = x_0.Size();

//  for (int time_step = 1; time_step <= n_time_steps; ++time_step)
//  {
//    const double cur_time = time_step * tstep;

//    const double r = Ricker(cur_time - tstep);

//    Vector y = x_1; y *= 2.0; y -= x_2;        // y = 2*x_1 - x_2
//    Vector z0; z0.SetSize(N); M.Mult(y, z0);   // z0 = M * (2*x_1 - x_2)
//    Vector z1; z1.SetSize(N); S.Mult(x_1, z1); // z1 = S * x_1
//    Vector z2 = *b; z2 *= r;                   // z2 = r * b
//    y = z1; y -= z2; y *= tstep*tstep;         // y = dt^2 * (S*x_1 - r*b)
//    Vector RHS = z0; RHS -= y;                 // RHS = M*(2*x_1-x_2) - dt^2*(S*x_1-r*b)

//    // Solve the SLAE
//    const int print_iter = 0;
//    PCG(M, prec, RHS, x_0, print_iter, 200, 1e-12, 0.0);

//    // Compute and print the L^2 norm of the error
//    if (time_step % tenth == 0)
//      cout << "step " << time_step << " / " << n_time_steps
//           << " ||solution||_{L^2} = " << x_0.Norml2() << endl;

//    Vector nodal_values_x, nodal_values_y;
//    x_0.GetNodalValues(nodal_values_x, 1);
//    x_0.GetNodalValues(nodal_values_y, 2);

//    if (visualization && (time_step % STEP_SNAP == 0))
//      save_vts(time_step, "displ", sx, sy, nx, ny, nodal_values_x, nodal_values_y);

//    for (int i = 0; i < nx+1; ++i)
//    {
//      float val_x = nodal_values_x(rec_y_index*(nx+1)+i);
//      float val_y = nodal_values_y(rec_y_index*(nx+1)+i);
//      solout_x.write((char*)&val_x, sizeof(val_x));
//      solout_y.write((char*)&val_y, sizeof(val_y));
//    }

//    x_2 = x_1;
//    x_1 = x_0;
//  }

//  cout << "Time loop is over" << endl;

//  delete b;
//  delete m;
//  delete s;
//  delete fespace;
//  delete fec;
//}




void test_polynomials()
{
  /*
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

  */
}
