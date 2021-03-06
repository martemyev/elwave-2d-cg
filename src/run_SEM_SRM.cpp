#include "elastic_wave2D.hpp"
#include "GLL_quadrature.hpp"
#include "parameters.hpp"
#include "receivers.hpp"
#include "source.hpp"

#include <fstream>

using namespace std;
using namespace mfem;

//#define OUTPUT_MASS_MATRIX

double mass_damp_weight(const mfem::Vector& point, const Parameters& param);
double stif_damp_weight(const mfem::Vector& point, const Parameters& param);

void output_coef_values(const Parameters& param,
                        RhoCoefficient &rho_coef,
                        RhoFuncCoefficient &rho_damp_coef,
                        LambdaFuncCoefficient &lambda_coef,
                        MuFuncCoefficient &mu_coef,
                        Mesh &mesh);



void ElasticWave2D::run_SEM_SRM()
{
  bool generate_edges = 1;
  Mesh mesh(param.grid.nx, param.grid.ny, Element::QUADRILATERAL,
            generate_edges, param.grid.sx, param.grid.sy);
  const int dim = mesh.Dimension();
  const int n_elements = param.grid.nx*param.grid.ny;
  MFEM_VERIFY(n_elements == mesh.GetNE(), "Unexpected number of mesh elements");

  FiniteElementCollection *fec = new H1_FECollection(param.order, dim);
  FiniteElementSpace fespace(&mesh, fec, dim); //, Ordering::byVDIM);
  cout << "Number of unknowns: " << fespace.GetVSize() << endl;

  RhoCoefficient        rho_coef(param);
  LambdaFuncCoefficient lambda_coef  (stif_damp_weight, param);
  MuFuncCoefficient     mu_coef      (stif_damp_weight, param);
  RhoFuncCoefficient    rho_damp_coef(mass_damp_weight, param);

  output_coef_values(param, rho_coef, rho_damp_coef, lambda_coef, mu_coef, mesh);

  IntegrationRule segment_GLL;
  create_segment_GLL_rule(param.order, segment_GLL);
  IntegrationRule quad_GLL(segment_GLL, segment_GLL);

  ElasticityIntegrator *elast_int = new ElasticityIntegrator(lambda_coef, mu_coef);
  elast_int->SetIntRule(&quad_GLL);
  BilinearForm stif(&fespace);
  stif.AddDomainIntegrator(elast_int);
  stif.Assemble();
  stif.Finalize();
  const SparseMatrix& S = stif.SpMat();

  VectorMassIntegrator *mass_int = new VectorMassIntegrator(rho_coef);
  mass_int->SetIntRule(&quad_GLL);
  BilinearForm mass(&fespace);
  mass.AddDomainIntegrator(mass_int);
  mass.Assemble();
  mass.Finalize();
  const SparseMatrix& M = mass.SpMat();

#if defined(OUTPUT_MASS_MATRIX)
  {
    ofstream mout("mass_mat.dat");
    mass.PrintMatlab(mout);
    cout << "M.nnz = " << M.NumNonZeroElems() << endl;
  }
#endif

  VectorMassIntegrator *damp_int = new VectorMassIntegrator(rho_damp_coef);
  damp_int->SetIntRule(&quad_GLL);
  BilinearForm dampM(&fespace);
  dampM.AddDomainIntegrator(damp_int);
  dampM.Assemble();
  dampM.Finalize();
  SparseMatrix& D = dampM.SpMat();
  double omega = 2.0*M_PI*param.source.frequency; // angular frequency
  D *= 0.5*param.dt*omega;

  LinearForm b(&fespace);
  if (param.source.plane_wave)
  {
    PlaneWaveSource plane_wave_source(dim, param);
    VectorDomainLFIntegrator *plane_wave_int =
        new VectorDomainLFIntegrator(plane_wave_source);
    plane_wave_int->SetIntRule(&quad_GLL);
    b.AddDomainIntegrator(plane_wave_int);
    b.Assemble();
  }
  else
  {
    if (!strcmp(param.source.type, "pointforce"))
    {
      VectorPointForce vector_point_force(dim, param);
      VectorDomainLFIntegrator *point_force_int =
          new VectorDomainLFIntegrator(vector_point_force);
      point_force_int->SetIntRule(&quad_GLL);
      b.AddDomainIntegrator(point_force_int);
      b.Assemble();
    }
    else if (!strcmp(param.source.type, "momenttensor"))
    {
      MomentTensorSource momemt_tensor_source(dim, param);
      VectorDomainLFIntegrator *moment_tensor_int =
          new VectorDomainLFIntegrator(momemt_tensor_source);
      moment_tensor_int->SetIntRule(&quad_GLL);
      b.AddDomainIntegrator(moment_tensor_int);
      b.Assemble();
    }
    else MFEM_ABORT("Unknown source type: " + string(param.source.type));
  }
  cout << "||b||_L2 = " << b.Norml2() << endl;

  Vector diagM; M.GetDiag(diagM); // mass matrix is diagonal
  Vector diagD; D.GetDiag(diagD); // damping matrix is diagonal

  const string method_name = "SEM_";

  int n_rec_sets = param.sets_of_receivers.size();
  ofstream *seisU = new ofstream[N_ELAST_COMPONENTS*n_rec_sets]; // for displacement
  ofstream *seisV = new ofstream[N_ELAST_COMPONENTS*n_rec_sets]; // for velocity
  for (int r = 0; r < n_rec_sets; ++r)
  {
    const string desc = param.sets_of_receivers[r]->description();
#if defined(MFEM_DEBUG)
    cout << desc << "\n";
    param.sets_of_receivers[r]->print_receivers(mesh);
#endif
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

  GridFunction u_0(&fespace); // displacement
  GridFunction u_1(&fespace);
  GridFunction u_2(&fespace);
  GridFunction v_1(&fespace); // velocity
  u_0 = 0.0;
  u_1 = 0.0;
  u_2 = 0.0;

  const int n_time_steps = ceil(param.T / param.dt);
  const int tenth = 0.1 * n_time_steps;

  const string snapshot_filebase = method_name + param.extra_string;
  const int N = u_0.Size();

  cout << "N time steps = " << n_time_steps
       << "\nTime loop..." << endl;

  // choose the time function for the source
  double (*temp_function)(const SourceParameters&, double) = nullptr;
  if (!strcmp(param.source.time_function, "auto"))
  {
    if (!strcmp(param.source.type, "pointforce"))
      temp_function = &RickerWavelet;
    else if (!strcmp(param.source.type, "momenttensor"))
      temp_function = &GaussFirstDerivative;
    else MFEM_ABORT("Unknown source type: " + string(param.source.type));
  }
  else if (!strcmp(param.source.time_function, "ricker"))
    temp_function = &RickerWavelet;
  else if (!strcmp(param.source.time_function, "firstgauss"))
    temp_function = &GaussFirstDerivative;
  else if (!strcmp(param.source.time_function, "gauss"))
    temp_function = &GaussFunction;
  else MFEM_ABORT("Unknown source time function: " +
                  string(param.source.time_function));

  for (int time_step = 1; time_step <= n_time_steps; ++time_step)
  {
    const double cur_time = time_step * param.dt;
    // the value of the time-dependent part of the source
    const double time_val = (*temp_function)(param.source, cur_time - param.dt);

    Vector y = u_1; y *= 2.0; y -= u_2;        // y = 2*u_1 - u_2

    Vector z0; z0.SetSize(N);                  // z0 = M * (2*u_1 - u_2)
    for (int i = 0; i < N; ++i) z0[i] = diagM[i] * y[i];

    Vector z1; z1.SetSize(N); S.Mult(u_1, z1); // z1 = S * u_1
    Vector z2 = b; z2 *= time_val;             // z2 = timeval*source

    // y = dt^2 * (S*u_1 - timeval*source), where it can be
    // y = dt^2 * (S*u_1 - ricker*pointforce) OR
    // y = dt^2 * (S*u_1 - ricker*planewave)  OR
    // y = dt^2 * (S*u_1 - gaussfirstderivative*momenttensor)
    y = z1; y -= z2; y *= param.dt*param.dt;

    // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-timeval*source)
    Vector RHS = z0; RHS -= y;

    for (int i = 0; i < N; ++i) y[i] = diagD[i] * u_2[i]; // y = D * u_2

    // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-timeval*source) + D*u_2
    RHS += y;

    // (M+D)*x_0 = M*(2*x_1-x_2) - dt^2*(S*x_1-r*b) + D*x_2
    for (int i = 0; i < N; ++i) u_0[i] = RHS[i] / (diagM[i]+diagD[i]);

    // velocity: v = du/dt, we use the central difference here
    v_1  = u_0;
    v_1 -= u_2;
    v_1 /= 2.0*param.dt;

    // Compute and print the L^2 norm of the error
    if (time_step % tenth == 0)
      cout << "step " << time_step << " / " << n_time_steps
           << " ||solution||_{L^2} = " << u_0.Norml2() << endl;

    if (time_step % param.step_snap == 0)
      output_snapshots(time_step, snapshot_filebase, param, u_0, v_1);

    output_seismograms(param, mesh, u_0, v_1, seisU, seisV);

    u_2 = u_1;
    u_1 = u_0;
  }

  delete[] seisV;
  delete[] seisU;

  cout << "Time loop is over" << endl;

  delete fec;
}



double mass_damp_weight(const Vector& point, const Parameters& param)
{
  const double x = point(0);
  const double y = point(1);
  const bool left   = (!strcmp(param.bc.left,   "abs") ? true : false);
  const bool right  = (!strcmp(param.bc.right,  "abs") ? true : false);
  const bool bottom = (!strcmp(param.bc.bottom, "abs") ? true : false);
  const bool top    = (!strcmp(param.bc.top,    "abs") ? true : false);

  const double X0 = 0.0;
  const double X1 = param.grid.sx;
  const double Y0 = 0.0;
  const double Y1 = param.grid.sy;
  const double layer = param.bc.damp_layer;
  const double power = param.bc.damp_power;

  // coef for the mass matrix in a damping region is computed
  // C_M = C_Mmax * x^p, where
  // p is typically 3,
  // x changes from 0 at the interface between damping and non-damping regions
  // to 1 at the boundary - the farthest damping layer
  // C_M in the non-damping region is 0

  double weight = 0.0;
  if (left && x - layer <= X0)
    weight += pow((X0-x+layer)/layer, power);
  else if (right && x + layer >= X1)
    weight += pow((x+layer-X1)/layer, power);

  if (bottom && y - layer <= Y0)
    weight += pow((Y0-y+layer)/layer, power);
  else if (top && y + layer >= Y1)
    weight += pow((y+layer-Y1)/layer, power);

  return weight;
}



double stif_damp_weight(const Vector& point, const Parameters& param)
{
  const double x = point(0);
  const double y = point(1);
  const bool left   = (!strcmp(param.bc.left,   "abs") ? true : false);
  const bool right  = (!strcmp(param.bc.right,  "abs") ? true : false);
  const bool bottom = (!strcmp(param.bc.bottom, "abs") ? true : false);
  const bool top    = (!strcmp(param.bc.top,    "abs") ? true : false);

  const double X0 = 0.0;
  const double X1 = param.grid.sx;
  const double Y0 = 0.0;
  const double Y1 = param.grid.sy;
  const double layer = param.bc.damp_layer;
  const double power = param.bc.damp_power+1;
  const double C0 = log(100.0);

  // coef for the stif matrix in a damping region is computed
  // C_K = exp(-C0*alpha(x)*k_inc*x), where
  // C0 = ln(100)
  // alpha(x) = a_Max * x^p
  // p is typically 3,
  // x changes from 0 to 1 (1 at the boundary - the farthest damping layer)
  // C_K in the non-damping region is 1

  double weight = 1.0;
  if (left && x - layer <= X0)
    weight *= exp(-C0*pow((X0-x+layer)/layer, power));
  else if (right && x + layer >= X1)
    weight *= exp(-C0*pow((x+layer-X1)/layer, power));

  if (bottom && y - layer <= Y0)
    weight *= exp(-C0*pow((Y0-y+layer)/layer, power));
  else if (top && y + layer >= Y1)
    weight *= exp(-C0*pow((y+layer-Y1)/layer, power));

  return weight;
}



void show_SRM_damp_weights(const Parameters& param)
{
  Vector mass_damp((param.grid.nx+1)*(param.grid.ny+1));
  Vector stif_damp((param.grid.nx+1)*(param.grid.ny+1));

  const double hx = param.grid.sx / param.grid.nx;
  const double hy = param.grid.sy / param.grid.ny;

  for (int iy = 0; iy < param.grid.ny+1; ++iy)
  {
    const double y = (iy == param.grid.ny ? param.grid.sy : iy*hy);
    for (int ix = 0; ix < param.grid.nx+1; ++ix)
    {
      const double x = (ix == param.grid.nx ? param.grid.sx : ix*hx);
      Vector point(2);
      point(0) = x;
      point(1) = y;
      const double md = mass_damp_weight(point, param);
      const double sd = stif_damp_weight(point, param);
      mass_damp(iy*(param.grid.nx+1)+ix) = md;
      stif_damp(iy*(param.grid.nx+1)+ix) = sd;
    }
  }

  string fname = "mass_damping_weights.vts";
  write_vts_scalar(fname, "mass_weights", param.grid.sx, param.grid.sy,
                   param.grid.nx, param.grid.ny, mass_damp);

  fname = "stif_damping_weights.vts";
  write_vts_scalar(fname, "stif_weights", param.grid.sx, param.grid.sy,
                   param.grid.nx, param.grid.ny, stif_damp);
}



void output_coef_values(const Parameters& param,
                        RhoCoefficient& rho_coef,
                        RhoFuncCoefficient& rho_damp_coef,
                        LambdaFuncCoefficient& lambda_coef,
                        MuFuncCoefficient& mu_coef,
                        Mesh& mesh)
{
  FiniteElementCollection *testfec = new H1_FECollection(1);
  FiniteElementSpace testfespace(&mesh, testfec);
  GridFunction test(&testfespace);
  Vector test_nodal;
  string fname;
  string extra = param.extra_string;

  test.ProjectCoefficient(rho_coef);
  test.GetNodalValues(test_nodal);
  fname = "rho_values_" + extra + ".vts";
  write_vts_scalar(fname, "rho", param.grid.sx, param.grid.sy,
                   param.grid.nx, param.grid.ny, test_nodal);

  test.ProjectCoefficient(rho_damp_coef);
  test.GetNodalValues(test_nodal);
  fname = "rho_damp_values_" + extra + ".vts";
  write_vts_scalar(fname, "rho_damp", param.grid.sx, param.grid.sy,
                   param.grid.nx, param.grid.ny, test_nodal);

  test.ProjectCoefficient(lambda_coef);
  test.GetNodalValues(test_nodal);
  fname = "lambda_values_" + extra + ".vts";
  write_vts_scalar(fname, "lambda", param.grid.sx, param.grid.sy,
                   param.grid.nx, param.grid.ny, test_nodal);

  test.ProjectCoefficient(mu_coef);
  test.GetNodalValues(test_nodal);
  fname = "mu_values_" + extra + ".vts";
  write_vts_scalar(fname, "mu", param.grid.sx, param.grid.sy,
                   param.grid.nx, param.grid.ny, test_nodal);

  delete testfec;
}

