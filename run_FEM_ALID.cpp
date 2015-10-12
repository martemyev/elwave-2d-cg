#include "elastic_wave2D.hpp"
#include "parameters.hpp"
#include "receivers.hpp"

#include <fstream>

using namespace std;
using namespace mfem;


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

void ElasticWave2D::run_FEM_ALID()
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

  rho_coef    = new CWConstCoefficient(param.rho_array, 0);
  lambda_coef = new CWConstCoefficient(lambda_array);
  mu_coef     = new CWConstCoefficient(mu_array);

  rho_damp_coef    = new CWFunctionCoefficient(mass_damp_weight, param, param.rho_array, 0);
  lambda_damp_coef = new CWFunctionCoefficient(mass_damp_weight, param, lambda_array, 0);
  mu_damp_coef     = new CWFunctionCoefficient(mass_damp_weight, param, mu_array, 0);

  elast_int = new ElasticityIntegrator(*lambda_coef, *mu_coef);
  stif = new BilinearForm(fespace);
  stif->AddDomainIntegrator(elast_int);

  mass_int = new VectorMassIntegrator(*rho_coef);
  mass = new BilinearForm(fespace);
  mass->AddDomainIntegrator(mass_int);

  double alpha1, alpha2;
  get_damp_alphas(param.source.frequency, alpha1, alpha2);

  dampS = new BilinearForm(fespace);
  dampS->AddDomainIntegrator(new ElasticityIntegrator(*lambda_damp_coef, *mu_damp_coef));
  dampS->Assemble();
  dampS->Finalize();
  SparseMatrix& D = dampS->SpMat();
  D *= alpha2;

  dampM = new BilinearForm(fespace);
  dampM->AddDomainIntegrator(new VectorMassIntegrator(*rho_damp_coef));
  dampM->Assemble();
  dampM->Finalize();
  SparseMatrix& DM = dampM->SpMat();
  DM *= alpha1;

  D += DM;
  D *= 0.5*param.dt;

  vector_point_force = new VectorPointForce(dim, param.source);
  point_force_int = new VectorDomainLFIntegrator(*vector_point_force);

  momemt_tensor_source = new MomentTensorSource(dim, param.source);
  moment_tensor_int = new VectorDomainLFIntegrator(*momemt_tensor_source);

  b = new LinearForm(fespace);
  b->AddDomainIntegrator(point_force_int);
  b->AddDomainIntegrator(moment_tensor_int);

  stif->Assemble();
  stif->Finalize();
  const SparseMatrix& S = stif->SpMat();

  mass->Assemble();
  mass->Finalize();
  SparseMatrix& M = mass->SpMat();

//  ofstream mout("mass_mat.dat");
//  mass->PrintMatlab(mout);
  cout << "M.nnz = " << M.NumNonZeroElems() << endl;

  SparseMatrix *Sys = nullptr;
  GSSmoother *prec = nullptr;

  const SparseMatrix& CopyFrom = D;
  const int nnz = CopyFrom.NumNonZeroElems();
  const bool ownij  = false;
  const bool ownval = true;
  Sys = new SparseMatrix(CopyFrom.GetI(), CopyFrom.GetJ(), new double[nnz],
                         CopyFrom.Height(), CopyFrom.Width(), ownij, ownval,
                         CopyFrom.areColumnsSorted());
  (*Sys) = 0.0;
  (*Sys) += D;
  (*Sys) += M;
  prec = new GSSmoother(*Sys);

  b->Assemble();
  cout << "||b||_L2 = " << b->Norml2() << endl;

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

    Vector z0; z0.SetSize(N); M.Mult(y, z0);   // z0 = M * (2*u_1 - u_2)

    Vector z1; z1.SetSize(N); S.Mult(u_1, z1); // z1 = S * u_1

    Vector z2 = *b; z2 *= r;                   // z2 = r * b

    y = z1; y -= z2; y *= param.dt*param.dt;   // y = dt^2 * (S*u_1 - r*b)

    Vector RHS = z0; RHS -= y;                 // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b)

    D.Mult(u_2, y);                            // y = D * u_2
    RHS += y;                                  // RHS = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b) + D*u_2
    // (M+D)*u_0 = M*(2*u_1-u_2) - dt^2*(S*u_1-r*b) + D*u_2
    PCG(*Sys, *prec, RHS, u_0, 0, 200, 1e-12, 0.0);

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

  delete prec;
  delete Sys;
}
