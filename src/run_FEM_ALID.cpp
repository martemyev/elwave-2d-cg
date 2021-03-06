#include "elastic_wave2D.hpp"
#include "parameters.hpp"
#include "receivers.hpp"
#include "source.hpp"

#include <fstream>

using namespace std;
using namespace mfem;

void compute_damping_weights(const Parameters& param, double *damping_weights);
void get_damp_alphas(double source_frequency, double &alpha1, double &alpha2);



void ElasticWave2D::run_FEM_ALID()
{
  bool generate_edges = 1;
  Mesh mesh(param.grid.nx, param.grid.ny, Element::QUADRILATERAL,
            generate_edges, param.grid.sx, param.grid.sy);
  const int dim = mesh.Dimension();
  const int n_elements = param.grid.nx*param.grid.ny;
  MFEM_VERIFY(n_elements == mesh.GetNE(), "Unexpected number of mesh elements");

  FiniteElementCollection *fec = new H1_FECollection(param.order, dim);
  FiniteElementSpace fespace(&mesh, fec, dim, Ordering::byVDIM);
  cout << "Number of unknowns: " << fespace.GetVSize() << endl;

  double *lambda_array = new double[n_elements];
  double *mu_array     = new double[n_elements];

  double *damp_weights = new double[n_elements];
  compute_damping_weights(param, damp_weights);

  double *rho_w_array   = new double[n_elements];
  double *lambda_w_array= new double[n_elements];
  double *mu_w_array    = new double[n_elements];

  for (int i = 0; i < n_elements; ++i)
  {
    const double rho = param.media.rho_array[i];
    const double vp  = param.media.vp_array[i];
    const double vs  = param.media.vs_array[i];
    const double w   = damp_weights[i];

    lambda_array[i]  = rho*(vp*vp - 2.*vs*vs);
    mu_array[i]      = rho*vs*vs;

    rho_w_array[i]   = rho*w;
    lambda_w_array[i]= lambda_array[i]*w;
    mu_w_array[i]    = mu_array[i]*w;
  }

  const bool own_array = false;
  CWConstCoefficient rho_coef(param.media.rho_array, own_array);
  CWConstCoefficient lambda_coef(lambda_array);
  CWConstCoefficient mu_coef(mu_array);

  CWConstCoefficient rho_w_coef(rho_w_array);
  CWConstCoefficient lambda_w_coef(lambda_w_array);
  CWConstCoefficient mu_w_coef(mu_w_array);

  BilinearForm stif(&fespace);
  stif.AddDomainIntegrator(new ElasticityIntegrator(lambda_coef, mu_coef));
  stif.Assemble();
  stif.Finalize();
  const SparseMatrix& S = stif.SpMat();

  BilinearForm mass(&fespace);
  mass.AddDomainIntegrator(new VectorMassIntegrator(rho_coef));
  mass.Assemble();
  mass.Finalize();
  SparseMatrix& M = mass.SpMat();

//  ofstream mout("mass_mat.dat");
//  mass->PrintMatlab(mout);
  cout << "M.nnz = " << M.NumNonZeroElems() << endl;

  double alpha1, alpha2;
  get_damp_alphas(param.source.frequency, alpha1, alpha2);

  BilinearForm dampS(&fespace);
  dampS.AddDomainIntegrator(new ElasticityIntegrator(lambda_w_coef, mu_w_coef));
  dampS.Assemble();
  dampS.Finalize();
  SparseMatrix& D = dampS.SpMat();
  D *= alpha2;

  BilinearForm dampM(&fespace);
  dampM.AddDomainIntegrator(new VectorMassIntegrator(rho_w_coef));
  dampM.Assemble();
  dampM.Finalize();
  SparseMatrix& DM = dampM.SpMat();
  DM *= alpha1;

  D += DM;
  D *= 0.5*param.dt;

  const SparseMatrix& CopyFrom = D;
  const int nnz = CopyFrom.NumNonZeroElems();
  const bool ownij  = false;
  const bool ownval = true;
  SparseMatrix Sys(CopyFrom.GetI(), CopyFrom.GetJ(), new double[nnz],
                   CopyFrom.Height(), CopyFrom.Width(), ownij, ownval,
                   CopyFrom.areColumnsSorted());
  Sys = 0.0;
  Sys += D;
  Sys += M;
  GSSmoother prec(Sys);

  VectorPointForce vector_point_force(dim, param);
  MomentTensorSource momemt_tensor_source(dim, param);

  LinearForm b(&fespace);
  b.AddDomainIntegrator(new VectorDomainLFIntegrator(vector_point_force));
  b.AddDomainIntegrator(new VectorDomainLFIntegrator(momemt_tensor_source));
  b.Assemble();
  cout << "||b||_L2 = " << b.Norml2() << endl;

  const string method_name = "FEM_";

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

  for (int time_step = 1; time_step <= n_time_steps; ++time_step)
  {
    const double cur_time = time_step * param.dt;

    const double r = RickerWavelet(param.source, cur_time - param.dt);

    Vector y = u_1; y *= 2.0; y -= u_2;        // y = 2*u_1 - u_2
    Vector z0; z0.SetSize(N); M.Mult(y, z0);   // z0 = M * (2*u_1 - u_2)
    Vector z1; z1.SetSize(N); S.Mult(u_1, z1); // z1 = S * u_1
    Vector z2 = b; z2 *= r;                    // z2 = r * b
    y = z1; y -= z2; y *= param.dt*param.dt;   // y = dt^2 * (S*u_1 - r*b)
    Vector RHS = z0; RHS -= y;                 // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b)
    D.Mult(u_2, y);                            // y = D * u_2
    RHS += y;                                  // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b) + D*u_2

    // (M+D)*u_0 = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b) + D*u_2
    PCG(Sys, prec, RHS, u_0, 0, 200, 1e-12, 0.0);

    // velocity
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



void compute_damping_weights(const Parameters& param, double *damping_weights)
{
  const bool left   = (!strcmp(param.bc.left,   "abs") ? true : false);
  const bool right  = (!strcmp(param.bc.right,  "abs") ? true : false);
  const bool bottom = (!strcmp(param.bc.bottom, "abs") ? true : false);
  const bool top    = (!strcmp(param.bc.top,    "abs") ? true : false);

  const int nx = param.grid.nx;
  const int ny = param.grid.ny;
  const double X0 = 0.0;
  const double X1 = param.grid.sx;
  const double Y0 = 0.0;
  const double Y1 = param.grid.sy;
  const double layer = param.bc.damp_layer;
  const double power = param.bc.damp_power;

  const double hx = (X1 - X0) / nx;
  const double hy = (Y1 - Y0) / ny;

  for (int ely = 0; ely < ny; ++ely)
  {
    const double y = Y0 + (ely+0.5)*hy; // center of a cell
    for (int elx = 0; elx < nx; ++elx)
    {
      const double x = X0 + (elx+0.5)*hx; // center of a cell

      double weight = 1e-12;

      if (left && x - layer < X0)
        weight += pow((X0-x+layer)/layer, power);
      else if (right && x + layer > X1)
        weight += pow((x+layer-X1)/layer, power);

      if (bottom && y - layer < Y0)
        weight += pow((Y0-y+layer)/layer, power);
      else if (top && y + layer > Y1)
        weight += pow((y+layer-Y1)/layer, power);

      const int el = ely*nx + elx;
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
