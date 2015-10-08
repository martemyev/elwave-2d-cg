#include "elastic_wave2D.hpp"
#include "GLL_quadrature.hpp"
#include "parameters.hpp"
#include "receivers.hpp"
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
  , damp(nullptr)
  , rho_coef(nullptr)
  , lambda_coef(nullptr)
  , mu_coef(nullptr)
  , rho_w_coef(nullptr)
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

  delete rho_w_coef;
  delete mu_coef;
  delete lambda_coef;
  delete rho_coef;

  delete momemt_tensor_source;
  delete vector_point_force;

  delete damp;
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
    damp_int->SetIntRule(quad_GLL);
    point_force_int->SetIntRule(quad_GLL);
    moment_tensor_int->SetIntRule(quad_GLL);
  }

  online_stage();

  delete quad_GLL;
  delete segment_GLL;
}

void compute_mass_damping_weights(int nx, int ny, double X0, double X1,
                                  double Y0, double Y1, double damping_layer,
                                  bool left, bool right, bool bottom, bool top,
                                  double *damping_weights)
{
  const double p = 1.2;

  const double hx = (X1 - X0) / nx;
  const double hy = (Y1 - Y0) / ny;

  for (int ely = 0; ely < ny; ++ely)
  {
    const double y = Y0 + (ely+0.5)*hy; // center of a cell
    for (int elx = 0; elx < nx; ++elx)
    {
      const double x = X0 + (elx+0.5)*hx; // center of a cell

      double weight = 0.0;

      if (left && x - damping_layer < X0)
        weight += pow((X0-(x-damping_layer))/damping_layer, p);
      else if (right && x + damping_layer > X1)
        weight += pow((x+damping_layer-X1)/damping_layer, p);

      if (bottom && y - damping_layer < Y0)
        weight += pow((Y0-(y-damping_layer))/damping_layer, p);
      else if (top && y + damping_layer > Y1)
        weight += pow((y+damping_layer-Y1)/damping_layer, p);

      const int el = ely*nx + elx;
      damping_weights[el] = weight;
    }
  }
}

void compute_stif_damping_weights(int nx, int ny, double X0, double X1,
                                  double Y0, double Y1, double damping_layer,
                                  bool left, bool right, bool bottom, bool top,
                                  double *damping_weights)
{
  const double p = 0.2;

  const double hx = (X1 - X0) / nx;
  const double hy = (Y1 - Y0) / ny;

  for (int ely = 0; ely < ny; ++ely)
  {
    const double y = Y0 + (ely+0.5)*hy; // center of a cell
    for (int elx = 0; elx < nx; ++elx)
    {
      const double x = X0 + (elx+0.5)*hx; // center of a cell

      double weight = 1.0;

      if (left && x - damping_layer < X0)
        weight *= exp(-log(100) * pow((X0-(x-damping_layer))/damping_layer, p+1));
      else if (right && x + damping_layer > X1)
        weight *= exp(-log(100) * pow((x+damping_layer-X1)/damping_layer, p+1));

      if (bottom && y - damping_layer < Y0)
        weight *= exp(-log(100) * pow((Y0-(y-damping_layer))/damping_layer, p+1));
      else if (top && y + damping_layer > Y1)
        weight *= exp(-log(100) * pow((y+damping_layer-Y1)/damping_layer, p+1));

      const int el = ely*nx + elx;
      damping_weights[el] = weight;
    }
  }
}

void get_damp_alpha(double source_frequency, double &alpha)
{
  // This is obtained by Shubin Fu, PhD student, Texas A&M
  const double q = 0.4;
  alpha = 2.0 * M_PI * source_frequency * q;
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
  double *mass_damp_weights = new double[n_elements];
  bool top_abs = (param.topsurf == 1 ? false : true);
  compute_mass_damping_weights(param.nx, param.ny, 0, param.sx, 0, param.sy,
                               param.damp_layer, 1, 1, 1, top_abs,
                               mass_damp_weights);
  double *stif_damp_weights = new double[n_elements];
  compute_stif_damping_weights(param.nx, param.ny, 0, param.sx, 0, param.sy,
                               param.damp_layer, 1, 1, 1, top_abs,
                               stif_damp_weights);
  double *rho_w_array = new double[n_elements];
  for (int i = 0; i < n_elements; ++i)
  {
    const double rho = param.rho_array[i];
    const double vp  = param.vp_array[i];
    const double vs  = param.vs_array[i];

    lambda_array[i]  = rho*(vp*vp - 2.*vs*vs);
    mu_array[i]      = rho*vs*vs;

    rho_w_array[i]   = rho*mass_damp_weights[i];

    lambda_array[i] *= stif_damp_weights[i];
    mu_array[i]     *= stif_damp_weights[i];
  }

  write_binary("stif_damp_weights.bin", n_elements, stif_damp_weights);
  write_binary("mass_damp_weights.bin", n_elements, mass_damp_weights);

  delete[] stif_damp_weights;
  delete[] mass_damp_weights;

  rho_coef = new CWConstCoefficient(param.rho_array, 0);
  lambda_coef = new CWConstCoefficient(lambda_array);
  mu_coef = new CWConstCoefficient(mu_array);
  rho_w_coef = new CWConstCoefficient(rho_w_array);

  elast_int = new ElasticityIntegrator(*lambda_coef, *mu_coef);
  stif = new BilinearForm(fespace);
  stif->AddDomainIntegrator(elast_int);

  mass_int = new VectorMassIntegrator(*rho_coef);
  mass = new BilinearForm(fespace);
  mass->AddDomainIntegrator(mass_int);

  damp_int = new VectorMassIntegrator(*rho_w_coef);
  damp = new BilinearForm(fespace);
  damp->AddDomainIntegrator(damp_int);

  vector_point_force = new VectorPointForce(dim, param.source);
  point_force_int = new VectorDomainLFIntegrator(*vector_point_force);

  momemt_tensor_source = new MomentTensorSource(dim, param.source);
  moment_tensor_int = new VectorDomainLFIntegrator(*momemt_tensor_source);

  b = new LinearForm(fespace);
  b->AddDomainIntegrator(point_force_int);
  b->AddDomainIntegrator(moment_tensor_int);
}

Vector compute_solution_at_points(const vector<Vertex>& points,
                                  const vector<int>& cells_containing_points,
                                  const GridFunction& U)
{
  MFEM_ASSERT(points.size() == cells_containing_points.size(), "Sizes mismatch");
  Vector U_at_points(2*points.size());
  Vector values(2);
  IntegrationPoint ip;
  for (size_t p = 0; p < points.size(); ++p)
  {
    ip.x = points[p](0);
    ip.y = points[p](1);
    ip.z = points[p](2);
    U.GetVectorValue(cells_containing_points[p], ip, values);
    MFEM_ASSERT(values.Size() == 2, "Unexpected vector size");
    U_at_points(2*p+0) = values(0);
    U_at_points(2*p+1) = values(1);
  }
  return U_at_points;
}

void ElasticWave2D::online_stage()
{
  stif->Assemble();
  stif->Finalize();
  const SparseMatrix& S = stif->SpMat();

  mass->Assemble();
  mass->Finalize();
  SparseMatrix& M = mass->SpMat();

  damp->Assemble();
  damp->Finalize();
  SparseMatrix& D = damp->SpMat();
  double alpha;
  get_damp_alpha(param.source.frequency, alpha);
  D *= alpha;
  D *= 0.5*param.dt;

//  ofstream mout("mass_mat.dat");
//  mass->PrintMatlab(mout);
  cout << "M.nnz = " << M.NumNonZeroElems() << endl;

  SparseMatrix *Sys = nullptr;
  GSSmoother *prec = nullptr;
  Vector *diagM = nullptr;
  Vector *diagD = nullptr;
  if (param.method == 0) // FEM
  {
    const SparseMatrix& CopyFrom = M;
    const int nnz = CopyFrom.NumNonZeroElems();
    const bool ownij  = false;
    const bool ownval = true;
    Sys = new SparseMatrix(CopyFrom.GetI(), CopyFrom.GetJ(), new double[nnz],
                           CopyFrom.Height(), CopyFrom.Width(), ownij, ownval,
                           CopyFrom.areColumnsSorted());
    (*Sys) = 0.0;
    (*Sys) += M;
    (*Sys) += D;
    prec = new GSSmoother(*Sys);
  }
  else if (param.method == 1) // SEM
  {
    diagM = new Vector;
    M.GetDiag(*diagM); // mass matrix is diagonal
    diagD = new Vector;
    D.GetDiag(*diagD); // damping matrix is diagonal
  }
  else
    mfem_error("Unknown method to be used");

  b->Assemble();
  cout << "||b||_L2 = " << b->Norml2() << endl;

  const string method_name = (param.method == 0 ? "FEM_" : "SEM_");

  int n_rec_sets = param.sets_of_receivers.size();
  ofstream *seisU = new ofstream[N_ELAST_COMPONENTS*n_rec_sets]; // for displacement
  ofstream *seisV = new ofstream[N_ELAST_COMPONENTS*n_rec_sets]; // for velocity
  for (int r = 0; r < n_rec_sets; ++r)
  {
    const string desc = param.sets_of_receivers[r]->description();
    for (int c = 0; c < N_ELAST_COMPONENTS; ++c)
    {
      string seismofile = method_name + param.extra_string + desc + "_u" + d2s(c) + ".bin";
      seisU[r*N_ELAST_COMPONENTS + c].open(seismofile.c_str(), ios::binary);
      MFEM_VERIFY(seisU[r*N_ELAST_COMPONENTS + c], "File '" + seismofile +
                  "' can't be opened");

      seismofile = method_name + param.extra_string + desc + "_v" + d2s(c) + ".bin";
      seisV[r*N_ELAST_COMPONENTS + c].open(seismofile.c_str(), ios::binary);
      MFEM_VERIFY(seisV[r*N_ELAST_COMPONENTS + c], "File '" + seismofile +
                  "' can't be opened");
    } // loop for components
  } // loop for sets of receivers

  GridFunction u_0(fespace); // displacement
  GridFunction u_1(fespace);
  GridFunction u_2(fespace);
  GridFunction v_1(fespace); // velocity
  u_0 = 0.0;
  u_1 = 0.0;
  u_2 = 0.0;

  Vector tmp;
  u_0.GetNodalValues(tmp, 1);
  cout << "nodal values size = " << tmp.Size() << endl;

  const int n_time_steps = ceil(param.T / param.dt);
  const int tenth = 0.1 * n_time_steps;

  const string snapshot_filebase = method_name + param.extra_string;
  const int N = u_0.Size();

  cout << "N time steps = " << n_time_steps
       << "\nTime loop..." << endl;

  int nbytes = 0;
  for (int time_step = 1; time_step <= n_time_steps; ++time_step)
  {
    const double cur_time = time_step * param.dt;

    const double r = param.source.Ricker(cur_time - param.dt);

    Vector y = u_1; y *= 2.0; y -= u_2;        // y = 2*u_1 - u_2

    Vector z0; z0.SetSize(N);                  // z0 = M * (2*u_1 - u_2)
    if (param.method == 0) // FEM
      M.Mult(y, z0);
    else if (param.method == 1) // SEM
      for (int i = 0; i < N; ++i) z0[i] = (*diagM)[i] * y[i];

    Vector z1; z1.SetSize(N); S.Mult(u_1, z1); // z1 = S * u_1

    Vector z2 = *b; z2 *= r;                   // z2 = r * b

    y = z1; y -= z2; y *= param.dt*param.dt;   // y = dt^2 * (S*u_1 - r*b)

    Vector RHS = z0; RHS -= y;                 // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b)

    if (param.method == 0) // FEM
    {
      D.Mult(u_2, y);                          // y = D * u_2
      RHS += y;                                // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b) + D*u_2
      // (M+D)*u_0 = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b) + D*u_2
      PCG(*Sys, *prec, RHS, u_0, 0, 200, 1e-12, 0.0);
    }
    else if (param.method == 1) // SEM
    {
      for (int i = 0; i < N; ++i) y[i] = (*diagD)[i] * u_2[i]; // y = D * u_2
      RHS += y;                                                // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b) + D*u_2
      // (M+D)*x_0 = M*(2*x_1-x_2) - dt^2*(S*x_1-r*b) + D*x_2
      for (int i = 0; i < N; ++i) u_0[i] = RHS[i] / ((*diagM)[i]+(*diagD)[i]);
    }

    // velocity
    v_1  = u_0;
    v_1 -= u_2;
    v_1 /= 2.0*param.dt;

    // Compute and print the L^2 norm of the error
    if (time_step % tenth == 0)
      cout << "step " << time_step << " / " << n_time_steps
           << " ||solution||_{L^2} = " << u_0.Norml2() << endl;

    if (time_step % param.step_snap == 0)
    {
      Vector u_x, u_y, v_x, v_y;
      u_0.GetNodalValues(u_x, 1);
      u_0.GetNodalValues(u_y, 2);
      v_1.GetNodalValues(v_x, 1);
      v_1.GetNodalValues(v_y, 2);

      string tstep = d2s(time_step,0,0,0,6);
      string fname = snapshot_filebase + "_U_t" + tstep + ".vts";
      write_vts(fname, "U", param.sx, param.sy, param.nx, param.ny, u_x, u_y);
      fname = snapshot_filebase + "_V_t" + tstep + ".vts";
      write_vts(fname, "V", param.sx, param.sy, param.nx, param.ny, v_x, v_y);
      fname = snapshot_filebase + "_Ux_t" + tstep + ".bin";
      write_binary(fname.c_str(), u_x.Size(), u_x);
      fname = snapshot_filebase + "_Uy_t" + tstep + ".bin";
      write_binary(fname.c_str(), u_y.Size(), u_y);
      fname = snapshot_filebase + "_Vx_t" + tstep + ".bin";
      write_binary(fname.c_str(), v_x.Size(), v_x);
      fname = snapshot_filebase + "_Vy_t" + tstep + ".bin";
      write_binary(fname.c_str(), v_y.Size(), v_y);
    }

    // for each set of receivers
    for (int rec = 0; rec < n_rec_sets; ++rec)
    {
      const ReceiversSet *rec_set = param.sets_of_receivers[rec];
      const Vector U_0 = compute_solution_at_points(rec_set->get_receivers(),
                                                    rec_set->get_cells_containing_receivers(),
                                                    u_0);
      const Vector U_2 = compute_solution_at_points(rec_set->get_receivers(),
                                                    rec_set->get_cells_containing_receivers(),
                                                    u_2);

      MFEM_ASSERT(U_0.Size() == N_ELAST_COMPONENTS*rec_set->n_receivers(),
                  "Sizes mismatch");
      Vector V_1 = U_0;
      V_1 -= U_2;
      V_1 /= 2.0*param.dt; // central difference

      float val;
      int n_loc_bytes = 0;
      for (int i = 0; i < U_0.Size(); i += N_ELAST_COMPONENTS)
      {
        for (int j = 0; j < N_ELAST_COMPONENTS; ++j)
        {
          val = U_0(i+j);
          seisU[rec*N_ELAST_COMPONENTS + j].write(reinterpret_cast<char*>(&val), sizeof(val));
//          cout << sizeof(val);

          val = V_1(i+j);
          seisV[rec*N_ELAST_COMPONENTS + j].write(reinterpret_cast<char*>(&val), sizeof(val));
        }
        n_loc_bytes += sizeof(float);
        nbytes += sizeof(float);
      }

      cout << n_loc_bytes << endl;
    } // for each set of receivers

    u_2 = u_1;
    u_1 = u_0;
  }


  cout << "nbytes = " << nbytes <<endl;
  cout << "Time loop is over" << endl;

  delete diagD;
  delete diagM;
  delete prec;
  delete Sys;
}




