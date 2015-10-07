#include "elastic_wave2D.hpp"
#include "GLL_quadrature.hpp"
#include "parameters.hpp"
#include "utilities.hpp"

#include <fstream>

using namespace std;
using namespace mfem;

#define RECEIVER_Y 700

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
VectorPointForce::VectorPointForce(int dim, const Source& s)
  : VectorCoefficient(dim)
  , source(s)
{ }

void VectorPointForce::Eval(Vector &V, ElementTransformation &T,
                            const IntegrationPoint &ip)
{
  double x[3];
  Vector transip(x, 3);
  T.Transform(ip, transip);
  V.SetSize(vdim);
  source.PointForce(transip, V);
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
MomentTensorSource::MomentTensorSource(int dim, const Source& s)
  : VectorCoefficient(dim)
  , source(s)
{ }

void MomentTensorSource::Eval(Vector &V, ElementTransformation &T,
                              const IntegrationPoint &ip)
{
  double x[3];
  Vector transip(x, 3);
  T.Transform(ip, transip);
  V.SetSize(vdim);
  source.MomentTensorSource(transip, V);
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
ElasticWave2D::ElasticWave2D(const Parameters &_param)
  : param(_param)
  , mesh(nullptr)
  , fec(nullptr)
  , fespace(nullptr)
  , stif(nullptr)
  , mass(nullptr)
  , rho_coef(nullptr)
  , lambda_coef(nullptr)
  , mu_coef(nullptr)
  , elast_int(nullptr)
  , mass_int(nullptr)
  , vector_point_force(nullptr)
  , momemt_tensor_source(nullptr)
  , point_force_int(nullptr)
  , moment_tensor_int(nullptr)
  , b(nullptr)
{ }

ElasticWave2D::~ElasticWave2D()
{
  delete b;

  delete mu_coef;
  delete lambda_coef;
  delete rho_coef;

  delete momemt_tensor_source;
  delete vector_point_force;

  delete mass;
  delete stif;
  delete fespace;
  delete fec;
  delete mesh;
}

void ElasticWave2D::run()
{
  offline_stage();

  IntegrationRule *segment_GLL = nullptr;
  IntegrationRule *quad_GLL = nullptr;

  if (param.method == 1) // SEM
  {
    segment_GLL = new IntegrationRule;
    create_segment_GLL_rule(param.order, *segment_GLL);
    quad_GLL = new IntegrationRule(*segment_GLL, *segment_GLL);

    elast_int->SetIntRule(quad_GLL);
    mass_int->SetIntRule(quad_GLL);
    point_force_int->SetIntRule(quad_GLL);
    moment_tensor_int->SetIntRule(quad_GLL);
  }

  online_stage();

  delete quad_GLL;
  delete segment_GLL;
}

void ElasticWave2D::offline_stage()
{
  bool generate_edges = 1;
  mesh = new Mesh(param.nx, param.ny, Element::QUADRILATERAL, generate_edges,
                  param.sx, param.sy);
  const int dim = mesh->Dimension();
  MFEM_VERIFY(param.nx*param.ny == mesh->GetNE(), "Unexpected number of mesh "
              "elements");

  fec = new H1_FECollection(param.order, dim);
  fespace = new FiniteElementSpace(mesh, fec, dim, Ordering::byVDIM);
  cout << "Number of unknowns: " << fespace->GetVSize() << endl;

  const int n_elements = param.nx*param.ny;
  double *lambda_array = new double[n_elements];
  double *mu_array     = new double[n_elements];
  for (int i = 0; i < n_elements; ++i)
  {
    const double rho = param.rho_array[i];
    const double vp  = param.vp_array[i];
    const double vs  = param.vs_array[i];

    lambda_array[i]  = rho*(vp*vp - 2.*vs*vs);
    mu_array[i]      = rho*vs*vs;
  }
  rho_coef = new CWConstCoefficient(param.rho_array, 0);
  lambda_coef = new CWConstCoefficient(lambda_array);
  mu_coef = new CWConstCoefficient(mu_array);

  elast_int = new ElasticityIntegrator(*lambda_coef, *mu_coef);
  stif = new BilinearForm(fespace);
  stif->AddDomainIntegrator(elast_int);

  mass_int = new VectorMassIntegrator(*rho_coef);
  mass = new BilinearForm(fespace);
  mass->AddDomainIntegrator(mass_int);

  vector_point_force = new VectorPointForce(dim, param.source);
  point_force_int = new VectorDomainLFIntegrator(*vector_point_force);

  momemt_tensor_source = new MomentTensorSource(dim, param.source);
  moment_tensor_int = new VectorDomainLFIntegrator(*momemt_tensor_source);

  b = new LinearForm(fespace);
  b->AddDomainIntegrator(point_force_int);
  b->AddDomainIntegrator(moment_tensor_int);
}

void ElasticWave2D::online_stage()
{
  stif->Assemble();
  stif->Finalize();
  const SparseMatrix& S = stif->SpMat();

  mass->Assemble();
  mass->Finalize();
  const SparseMatrix& M = mass->SpMat();

//  ofstream mout("mass_mat.dat");
//  mass->PrintMatlab(mout);
  cout << "M.nnz = " << M.NumNonZeroElems() << endl;

  GSSmoother *prec = nullptr;
  Vector *diagM = nullptr;
  if (param.method == 0) // FEM
    prec = new GSSmoother(M);
  else if (param.method == 1) // SEM
  {
    diagM = new Vector;
    M.GetDiag(*diagM); // mass matrix is diagonal
  }
  else
    mfem_error("Unknown method to be used");

  b->Assemble();
  cout << "||b||_L2 = " << b->Norml2() << endl;

  const double hy = param.sy / param.ny;
  const int rec_y_index = RECEIVER_Y / hy;

  ofstream solout_x("seis0.bin", std::ios::binary);
  ofstream solout_y("seis1.bin", std::ios::binary);

  GridFunction x_0(fespace);
  GridFunction x_1(fespace);
  GridFunction x_2(fespace);
  x_0 = 0.0;
  x_1 = 0.0;
  x_2 = 0.0;

  Vector tmp;
  x_0.GetNodalValues(tmp, 1);
  cout << "nodal values size = " << tmp.Size() << endl;

  const int n_time_steps = ceil(param.T / param.dt);
  const int tenth = 0.1 * n_time_steps;

  const string method_name = (param.method == 0 ? "FEM" : "SEM");
  const string snapshot_filebase = method_name + "_o" + d2s(param.order);
  const int N = x_0.Size();

  cout << "N time steps = " << n_time_steps
       << "\nTime loop..." << endl;

  for (int time_step = 1; time_step <= n_time_steps; ++time_step)
  {
    const double cur_time = time_step * param.dt;

    const double r = param.source.Ricker(cur_time - param.dt);

    Vector y = x_1; y *= 2.0; y -= x_2;        // y = 2*x_1 - x_2

    Vector z0; z0.SetSize(N);                  // z0 = M * (2*x_1 - x_2)
    if (param.method == 0) // FEM
      M.Mult(y, z0);
    else if (param.method == 1) // SEM
      for (int i = 0; i < N; ++i) z0[i] = (*diagM)[i] * y[i];

    Vector z1; z1.SetSize(N); S.Mult(x_1, z1); // z1 = S * x_1

    Vector z2 = *b; z2 *= r;                   // z2 = r * b

    y = z1; y -= z2; y *= param.dt*param.dt;   // y = dt^2 * (S*x_1 - r*b)

    Vector RHS = z0; RHS -= y;                 // RHS = M*(2*x_1-x_2) - dt^2*(S*x_1-r*b)
    if (param.method == 0) // FEM
      PCG(M, *prec, RHS, x_0, 0, 200, 1e-12, 0.0);
    else if (param.method == 1) // SEM
      for (int i = 0; i < N; ++i) x_0[i] = RHS[i] / (*diagM)[i];

    // Compute and print the L^2 norm of the error
    if (time_step % tenth == 0)
      cout << "step " << time_step << " / " << n_time_steps
           << " ||solution||_{L^2} = " << x_0.Norml2() << endl;

    Vector nodal_values_x, nodal_values_y;
    x_0.GetNodalValues(nodal_values_x, 1);
    x_0.GetNodalValues(nodal_values_y, 2);

    if (time_step % param.step_snap == 0)
      save_vts(snapshot_filebase, time_step, "U", param.sx, param.sy, param.nx,
               param.ny, nodal_values_x, nodal_values_y);

    for (int i = 0; i < param.nx+1; ++i)
    {
      float val_x = nodal_values_x(rec_y_index*(param.nx+1)+i);
      float val_y = nodal_values_y(rec_y_index*(param.nx+1)+i);
      solout_x.write((char*)&val_x, sizeof(val_x));
      solout_y.write((char*)&val_y, sizeof(val_y));
    }

    x_2 = x_1;
    x_1 = x_0;
  }

  cout << "Time loop is over" << endl;

  delete diagM;
  delete prec;
}
