#include "mfem.hpp"
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
const double SOURCE_ANGLE  = 270 * M_PI/180.0;

const double RECEIVER_Y = 700;

const int STEP_SNAP = 500;

void Delta(const Vector&, Vector&);
void Gauss(const Vector&, Vector&);
double Ricker(double);

void save_vts(int time_step, const std::string &solname,
              double sx, double sy, int nx, int ny,
              const Vector& sol_x, const Vector& sol_y);


/**
 * Cell-wise constant coefficient
 */
class CWConstCoefficient : public Coefficient
{
public:
  CWConstCoefficient(double *array) : val_array(array) { }
  ~CWConstCoefficient() { delete[] val_array; }

  double Eval(ElementTransformation &T, const IntegrationPoint &ip)
  {
    return val_array[T.ElementNo];
  }

private:
  double *val_array;
};


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

void read_binary(const char *filename, int n_values, double *values);

void get_damp_alphas(double source_frequency, double &alpha1, double &alpha2);




int main(int argc, char *argv[])
{
   double sx = 1000.0, sy = 1000.0; // size of the computational domain, m
   int nx = 250, ny = 250; // number of elements in each direction
   double tend = 0.4; // simulation time, s
   double tstep = 1e-4; // time discretization step, s
   int order = 1; // order of the finite element basis functions for the approximation
   bool visualization = 1; // whether we output the solution for visualization
   const char *rhofile = "no-file", *vpfile = "no-file", *vsfile = "no-file";

   // Parse command-line options.
   OptionsParser args(argc, argv);
   args.AddOption(&sx, "-sx", "--sizex", "Size of domain in x-direction, m");
   args.AddOption(&sy, "-sy", "--sizey", "Size of domain in y-direction, m");
   args.AddOption(&nx, "-nx", "--numberx", "Number of elements in x-direction");
   args.AddOption(&ny, "-ny", "--numbery", "Number of elements in y-direction");
   args.AddOption(&tend, "-tend", "--time-end", "Simulation time, s");
   args.AddOption(&tstep, "-tstep", "--time-step", "Time step, s");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree)");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization", "Enable or disable GLVis visualization.");
   args.AddOption(&rhofile, "-rhofile", "--rhofile", "Density file");
   args.AddOption(&vpfile,  "-vpfile",  "--vpfile",  "P-wave velocity file");
   args.AddOption(&vsfile,  "-vsfile",  "--vsfile",  "S-wave velocity file");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   bool gen_edges = 1;
   Mesh *mesh = new Mesh(nx, ny, Element::QUADRILATERAL, gen_edges, sx, sy);
   int dim = mesh->Dimension();


   double *rho_array = new double[nx*ny];
   double *vp_array  = new double[nx*ny];
   double *vs_array  = new double[nx*ny];

//   double min_vp, max_vp, min_vs, max_vs;

   if (strcmp(rhofile, "no-file") == 0)
     for (int i = 0; i < nx*ny; ++i) rho_array[i] = RHO;
   else
     read_binary(rhofile, nx*ny, rho_array);

   if (strcmp(vpfile, "no-file") == 0)
     for (int i = 0; i < nx*ny; ++i) vp_array[i] = VP;
   else
     read_binary(vpfile,  nx*ny, vp_array);

   if (strcmp(vsfile, "no-file") == 0)
     for (int i = 0; i < nx*ny; ++i) vs_array[i] = VS;
   else
     read_binary(vsfile,  nx*ny, vs_array);

   const double abs_layer_width = 100.0;
   bool left(true), right(true), bottom(true), top(false);
   double *damping_weights = new double[nx*ny];
   compute_rayleigh_damping_weights(nx, ny, 0, sx, 0, sy, abs_layer_width,
                                    left, right, bottom, top, damping_weights);

   double alpha1, alpha2;
   get_damp_alphas(FREQUENCY, alpha1, alpha2);
   cout << "alpha1 = " << alpha1 << "\n"
        << "alpha2 = " << alpha2 << "\n" << endl;

   double *lambda_array   = new double[nx*ny];
   double *mu_array       = new double[nx*ny];
   double *w_rho_array    = new double[nx*ny];
   double *w_lambda_array = new double[nx*ny];
   double *w_mu_array     = new double[nx*ny];
   for (int i = 0; i < nx*ny; ++i)
   {
     const double rho = rho_array[i];
     const double vp  = vp_array[i];
     const double vs  = vs_array[i];
     const double w   = damping_weights[i];

     lambda_array[i]  = rho*(vp*vp - 2.*vs*vs);
     mu_array[i]      = rho*vs*vs;
     w_rho_array[i]   = w*rho;
     w_lambda_array[i]= w*lambda_array[i];
     w_mu_array[i]    = w*mu_array[i];
   }
   CWConstCoefficient rho_cw(rho_array);
   CWConstCoefficient lambda_cw(lambda_array);
   CWConstCoefficient mu_cw(mu_array);
   CWConstCoefficient w_rho(w_rho_array);
   CWConstCoefficient w_lambda(w_lambda_array);
   CWConstCoefficient w_mu(w_mu_array);

   FiniteElementCollection *fec;
   fec = new H1_FECollection(order, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec, dim, Ordering::byNODES);
   cout << "Number of unknowns: " << fespace->GetVSize() << endl;


   BilinearForm *s = new BilinearForm(fespace);
   s->AddDomainIntegrator(new ElasticityIntegrator(lambda_cw, mu_cw));
   s->Assemble();
   s->Finalize();
   const SparseMatrix& S = s->SpMat();

   BilinearForm *m = new BilinearForm(fespace);
   m->AddDomainIntegrator(new VectorMassIntegrator(rho_cw));
   m->Assemble();
   m->Finalize();
   SparseMatrix& M = m->SpMat();


   BilinearForm *ds = new BilinearForm(fespace);
   ds->AddDomainIntegrator(new ElasticityIntegrator(w_lambda, w_mu));
   ds->Assemble();
   ds->Finalize();
   SparseMatrix& D = ds->SpMat();
   D *= alpha2;

   BilinearForm *dm = new BilinearForm(fespace);
   dm->AddDomainIntegrator(new VectorMassIntegrator(w_rho));
   dm->Assemble();
   dm->Finalize();
   SparseMatrix& DM = dm->SpMat();
   DM *= alpha1;

   D += DM;
   D *= 0.5*tstep;

   cout << "D.nnz = " << D.NumNonZeroElems() << "\n"
        << "M.nnz = " << M.NumNonZeroElems() << endl;


   const int nnz = D.NumNonZeroElems();
   const bool ownij  = false;
   const bool ownval = true;
   SparseMatrix Sys(D.GetI(), D.GetJ(), new double[nnz], D.Height(), D.Width(),
                    ownij, ownval, D.areColumnsSorted());
   Sys = 0.0;
   Sys += D;
   Sys += M;


   // Define a simple symmetric Gauss-Seidel preconditioner.
   GSSmoother prec(Sys);

   GridFunction x_0(fespace);
   GridFunction x_1(fespace);
   GridFunction x_2(fespace);
   x_0 = 0.0;
   x_1 = 0.0;
   x_2 = 0.0;


   // Set up the linear form b(.) which corresponds to the RHS.
   LinearForm *b = new LinearForm(fespace);
   VectorFunctionCoefficient f(dim, Gauss);
   b->AddDomainIntegrator(new VectorDomainLFIntegrator(f));
   b->Assemble();

   cout << "||b||_L2 = " << b->Norml2() << endl;

   const double hy = sy / ny;
   const int rec_y_index = RECEIVER_Y / hy;

   ofstream solout_x("seis0.bin", std::ios::binary);
   ofstream solout_y("seis1.bin", std::ios::binary);

   const int n_time_steps = ceil(tend/tstep);
   const int tenth = 0.1 * n_time_steps;

   cout << "N time steps = " << n_time_steps
        << "\nTime loop..." << endl;

   const int N = x_0.Size();

   for (int time_step = 1; time_step <= n_time_steps; ++time_step)
   {
     const double cur_time = time_step * tstep;

     const double r = Ricker(cur_time - tstep);

     Vector y = x_1; y *= 2.0; y -= x_2;        // y = 2*x_1 - x_2
     Vector z0; z0.SetSize(N); M.Mult(y, z0);   // z0 = M * (2*x_1 - x_2)
     Vector z1; z1.SetSize(N); S.Mult(x_1, z1); // z1 = S * x_1
     Vector z2 = *b; z2 *= r;                   // z2 = r * b
     y = z1; y -= z2; y *= tstep*tstep;         // y = dt^2 * (S*x_1 - r*b)
     Vector RHS = z0; RHS -= y;                 // RHS = M*(2*x_1-x_2) - dt^2*(S*x_1-r*b)
     D.Mult(x_2, y);                            // y = D * x_2
     RHS += y;                                  // RHS = M*(2*x_1-x_2) - dt^2*(S*x_1-r*b) + D*x_2

     // Solve the SLAE
     const int print_iter = 0;
     PCG(Sys, prec, RHS, x_0, print_iter, 200, 1e-12, 0.0);

     // Compute and print the L^2 norm of the error
     if (time_step % tenth == 0)
       cout << "step " << time_step << " / " << n_time_steps
            << " ||solution||_{L^2} = " << x_0.Norml2() << endl;

     Vector nodal_values_x, nodal_values_y;
     x_0.GetNodalValues(nodal_values_x, 1);
     x_0.GetNodalValues(nodal_values_y, 2);

     if (visualization && (time_step % STEP_SNAP == 0))
     {
       save_vts(time_step, "displ", sx, sy, nx, ny, nodal_values_x, nodal_values_y);
     }

     for (int i = 0; i < nx+1; ++i)
     {
       float val_x = nodal_values_x(rec_y_index*(nx+1)+i);
       float val_y = nodal_values_y(rec_y_index*(nx+1)+i);
       solout_x.write((char*)&val_x, sizeof(val_x));
       solout_y.write((char*)&val_y, sizeof(val_y));
     }

     x_2 = x_1;
     x_1 = x_0;
   }

   cout << "Time loop is over" << endl;

   // Free the used memory.
   delete b;
   delete dm;
   delete ds;
   delete m;
   delete s;
   delete fespace;
   delete fec;
   delete mesh;
}



void Delta(const Vector& x, Vector& f)
{
  double center[] = {SOURCE_X, SOURCE_Y}; // source location
  const double tol = 1e-2;
  if (x.DistanceTo(center) < tol)
  {
    f(0) = cos(SOURCE_ANGLE);
    f(1) = sin(SOURCE_ANGLE);
  }
  else
    f = 0.0;
}

void Gauss(const Vector& x, Vector& f)
{
  const double G = exp(-((x(0)-SOURCE_X)*(x(0)-SOURCE_X) + (x(1)-SOURCE_Y)*(x(1)-SOURCE_Y)) / GAUSS_SUPPORT / GAUSS_SUPPORT);
  f(0) = G*cos(SOURCE_ANGLE);
  f(1) = G*sin(SOURCE_ANGLE);
}

double Ricker(double t)
{
  const double part  = M_PI*FREQUENCY*(t-1./FREQUENCY);
  const double part2 = part*part;
  return RICKER_SCALE * (1. - 2.*part2)*exp(-part2);
}

void save_vts(int time_step, const std::string &solname,
              double sx, double sy, int nx, int ny,
              const Vector& sol_x, const Vector& sol_y)
{
  ostringstream ostr;
  ostr << time_step;
  const string filename = "results_Dleft_t" + ostr.str() + ".vts";
  ofstream out(filename.c_str());
  if (!out) throw runtime_error("File can't be opened");

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
  out << "  <StructuredGrid WholeExtent=\"1 " << nx+1 << " 1 1 1 " << ny+1 << "\">\n";
  out << "    <Piece Extent=\"1 " << nx+1 << " 1 1 1 " << ny+1 << "\">\n";
  out << "      <PointData Vectors=\"" << solname << "\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"" << solname
      << "\" format=\"ascii\" NumberOfComponents=\"3\">\n";

  for (int i = 0; i < ny+1; ++i)
  {
    for (int j = 0; j < nx+1; ++j)
    {
      const int glob_vert_index = i*(nx+1) + j;
      out << sol_x(glob_vert_index) << " "
          << sol_y(glob_vert_index) << " 0.0 ";
    }
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </PointData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         "format=\"ascii\">\n";

  const double hx = sx / nx;
  const double hy = sy / ny;

  for (int i = 0; i < ny+1; ++i)
  {
    const double y = (i == ny ? sy : i*hy);
    for (int j = 0; j < nx+1; ++j)
    {
      const double x = (j == nx ? sx : j*hx);
      out << x << " " << y << " 0.0 ";
    }
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "    </Piece>\n";
  out << "  </StructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
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



void read_binary(const char *filename,
                 int n_values,
                 double *values)
{
  std::ifstream in(filename, std::ios::binary);
  MFEM_VERIFY(in, "File '" + string(filename) + "' can't be opened");

  in.seekg(0, in.end); // jump to the end of the file
  int length = in.tellg(); // total length of the file in bytes
  int size_value = length / n_values; // size (in bytes) of one value

  MFEM_VERIFY(length % n_values != 0, "The number of bytes in the file '" +
              string(filename) + "' is not divisible by the number of elements "
              + d2s(n_values));

  in.seekg(0, in.beg); // jump to the beginning of the file

  if (size_value == sizeof(double))
  {
    in.read((char*)values, n_values*size_value); // read all at once

    MFEM_VERIFY(n_values == static_cast<int>(in.gcount()), "The number of "
                "successfully read elements is different from the expected one");
  }
  else if (size_value == sizeof(float))
  {
    float val = 0;
    for (int i = 0; i < n_values; ++i)  // read element-by-element
    {
      in.read((char*)&val, size_value); // read a 'float' value
      values[i] = val;                  // convert it to a 'double' value
    }
  }
  else MFEM_VERIFY(0, "Unknown size of an element (" + d2s(size_value) + ") in "
                   "bytes. Expected one is either sizeof(float) = " +
                   d2s(sizeof(float)) + ", or sizeof(double) = " +
                   d2s(sizeof(double)));

  in.close();
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
