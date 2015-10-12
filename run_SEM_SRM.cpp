#include "elastic_wave2D.hpp"
#include "GLL_quadrature.hpp"
#include "parameters.hpp"
#include "receivers.hpp"

#include <fstream>

using namespace std;
using namespace mfem;



void ElasticWave2D::run_SEM_SRM()
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

  rho_coef      = new CWConstCoefficient(param.rho_array, 0);
  lambda_coef   = new CWFunctionCoefficient(stif_damp_weight, param, lambda_array);
  mu_coef       = new CWFunctionCoefficient(stif_damp_weight, param, mu_array);
  rho_damp_coef = new CWFunctionCoefficient(mass_damp_weight, param, param.rho_array, 0);

  elast_int = new ElasticityIntegrator(*lambda_coef, *mu_coef);
  stif = new BilinearForm(fespace);
  stif->AddDomainIntegrator(elast_int);

  mass_int = new VectorMassIntegrator(*rho_coef);
  mass = new BilinearForm(fespace);
  mass->AddDomainIntegrator(mass_int);

  damp_int = new VectorMassIntegrator(*rho_damp_coef);
  dampM = new BilinearForm(fespace);
  dampM->AddDomainIntegrator(damp_int);

  vector_point_force = new VectorPointForce(dim, param.source);
  point_force_int = new VectorDomainLFIntegrator(*vector_point_force);

  momemt_tensor_source = new MomentTensorSource(dim, param.source);
  moment_tensor_int = new VectorDomainLFIntegrator(*momemt_tensor_source);

  b = new LinearForm(fespace);
  b->AddDomainIntegrator(point_force_int);
  b->AddDomainIntegrator(moment_tensor_int);

  IntegrationRule *segment_GLL = nullptr;
  IntegrationRule *quad_GLL = nullptr;

  segment_GLL = new IntegrationRule;
  create_segment_GLL_rule(param.order, *segment_GLL);
  quad_GLL = new IntegrationRule(*segment_GLL, *segment_GLL);

  elast_int->SetIntRule(quad_GLL);
  mass_int->SetIntRule(quad_GLL);
  damp_int->SetIntRule(quad_GLL);
  point_force_int->SetIntRule(quad_GLL);
  moment_tensor_int->SetIntRule(quad_GLL);

  stif->Assemble();
  stif->Finalize();
  const SparseMatrix& S = stif->SpMat();

  mass->Assemble();
  mass->Finalize();
  SparseMatrix& M = mass->SpMat();

  dampM->Assemble();
  dampM->Finalize();
  SparseMatrix& D = dampM->SpMat();
  double omega = 2.0*M_PI*param.source.frequency; // angular frequency
  D *= 0.5*param.dt*omega;

//  ofstream mout("mass_mat.dat");
//  mass->PrintMatlab(mout);
  cout << "M.nnz = " << M.NumNonZeroElems() << endl;

  Vector *diagM = nullptr;
  Vector *diagD = nullptr;

  diagM = new Vector;
  M.GetDiag(*diagM); // mass matrix is diagonal
  diagD = new Vector;
  D.GetDiag(*diagD); // damping matrix is diagonal

  b->Assemble();
  cout << "||b||_L2 = " << b->Norml2() << endl;

  const string method_name = "SEM_";

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

  for (int time_step = 1; time_step <= n_time_steps; ++time_step)
  {
    const double cur_time = time_step * param.dt;

    const double r = param.source.Ricker(cur_time - param.dt);

    Vector y = u_1; y *= 2.0; y -= u_2;        // y = 2*u_1 - u_2

    Vector z0; z0.SetSize(N);                  // z0 = M * (2*u_1 - u_2)
    for (int i = 0; i < N; ++i) z0[i] = (*diagM)[i] * y[i];

    Vector z1; z1.SetSize(N); S.Mult(u_1, z1); // z1 = S * u_1

    Vector z2 = *b; z2 *= r;                   // z2 = r * b

    y = z1; y -= z2; y *= param.dt*param.dt;   // y = dt^2 * (S*u_1 - r*b)

    Vector RHS = z0; RHS -= y;                 // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b)

    for (int i = 0; i < N; ++i) y[i] = (*diagD)[i] * u_2[i]; // y = D * u_2
    RHS += y;                                                // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b) + D*u_2
    // (M+D)*x_0 = M*(2*x_1-x_2) - dt^2*(S*x_1-r*b) + D*x_2
    for (int i = 0; i < N; ++i) u_0[i] = RHS[i] / ((*diagM)[i]+(*diagD)[i]);

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
      write_vts_vector(fname, "U", param.sx, param.sy, param.nx, param.ny, u_x, u_y);
      fname = snapshot_filebase + "_V_t" + tstep + ".vts";
      write_vts_vector(fname, "V", param.sx, param.sy, param.nx, param.ny, v_x, v_y);
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
      for (int i = 0; i < U_0.Size(); i += N_ELAST_COMPONENTS)
      {
        for (int j = 0; j < N_ELAST_COMPONENTS; ++j)
        {
          val = U_0(i+j);
          seisU[rec*N_ELAST_COMPONENTS + j].write(reinterpret_cast<char*>(&val), sizeof(val));

          val = V_1(i+j);
          seisV[rec*N_ELAST_COMPONENTS + j].write(reinterpret_cast<char*>(&val), sizeof(val));
        }
      }
    } // for each set of receivers

    u_2 = u_1;
    u_1 = u_0;
  }

  delete[] seisV;
  delete[] seisU;

  cout << "Time loop is over" << endl;

  delete diagD;
  delete diagM;

  delete quad_GLL;
  delete segment_GLL;
}

