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
  CWConstCoefficient(double *array, bool own = 1)
    : val_array(array), own_array(own)
  { }

  ~CWConstCoefficient() { if (own_array) delete[] val_array; }

  double Eval(ElementTransformation &T, const IntegrationPoint &ip)
  {
    return val_array[T.ElementNo];
  }

private:
  double *val_array;
  bool own_array;
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
void get_minmax(double *array, int n_elements, double &min_val, double &max_val);

void get_damp_alphas(double source_frequency, double &alpha1, double &alpha2);

void solve_without_damp_layer(int nx, int ny, double *rho_array,
                              double *vp_array, double *vs_array, int order,
                              int dim, Mesh *mesh, double sx, double sy,
                              double tend, double tstep, bool visualization);




int main(int argc, char *argv[])
{
   double sx = 1000.0, sy = 1000.0; // size of the computational domain, m
   int nx = 250, ny = 250; // number of elements in each direction
   double tend = 0.4; // simulation time, s
   double tstep = 1e-4; // time discretization step, s
   int order = 1; // order of the finite element basis functions for the approximation
   bool visualization = 1; // whether we output the solution for visualization
   const char *rhofile = "no-file", *vpfile = "no-file", *vsfile = "no-file";
   const double abs_layer_width = 100.0;

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
   cout << "dim = " << dim << endl;
   const int n_elements = mesh->GetNE();
   MFEM_VERIFY(nx*ny == n_elements, "Unexpected number of mesh elements");

   double *rho_array = new double[n_elements];
   double *vp_array  = new double[n_elements];
   double *vs_array  = new double[n_elements];

   double min_rho, max_rho, min_vp, max_vp, min_vs, max_vs;

   if (strcmp(rhofile, "no-file") == 0)
   {
     for (int i = 0; i < n_elements; ++i) rho_array[i] = RHO;
     min_rho = max_rho = RHO;
   }
   else
   {
     read_binary(rhofile, n_elements, rho_array);
     get_minmax(rho_array, n_elements, min_rho, max_rho);
   }

   if (strcmp(vpfile, "no-file") == 0)
   {
     for (int i = 0; i < n_elements; ++i) vp_array[i] = VP;
     min_vp = max_vp = VP;
   }
   else
   {
     read_binary(vpfile, n_elements, vp_array);
     get_minmax(vp_array, n_elements, min_vp, max_vp);
   }

   if (strcmp(vsfile, "no-file") == 0)
   {
     for (int i = 0; i < n_elements; ++i) vs_array[i] = VS;
     min_vs = max_vs = VS;
   }
   else
   {
     read_binary(vsfile, n_elements, vs_array);
     get_minmax(vs_array, n_elements, min_vs, max_vs);
   }

   const double min_wavelength = min(min_vp, min_vs) / (2.0*FREQUENCY);
   cout << "min wavelength = " << min_wavelength << endl;

   if (abs_layer_width < 2.5*min_wavelength)
     mfem_warning("damping layer for absorbing bc should be about 3*wavelength");


//   solve_with_damp_layer();
   solve_without_damp_layer(nx, ny, rho_array,vp_array, vs_array, order, dim,
                            mesh, sx, sy, tend, tstep, visualization);

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



void get_minmax(double *array, int n_elements, double &min_val, double &max_val)
{
  min_val = max_val = array[0];
  for (int i = 1; i < n_elements; ++i)
  {
    min_val = min(min_val, array[i]);
    max_val = max(max_val, array[i]);
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





void solve_without_damp_layer(int nx, int ny, double *rho_array,
                              double *vp_array, double *vs_array, int order,
                              int dim, Mesh *mesh, double sx, double sy,
                              double tend, double tstep, bool visualization)
{
  const int n_elements = nx*ny;
  double *lambda_array   = new double[n_elements];
  double *mu_array       = new double[n_elements];
  for (int i = 0; i < n_elements; ++i)
  {
    const double rho = rho_array[i];
    const double vp  = vp_array[i];
    const double vs  = vs_array[i];

    lambda_array[i]  = rho*(vp*vp - 2.*vs*vs);
    mu_array[i]      = rho*vs*vs;
  }
  CWConstCoefficient rho_cw(rho_array);
  CWConstCoefficient lambda_cw(lambda_array);
  CWConstCoefficient mu_cw(mu_array);

  FiniteElementCollection *fec;
  fec = new H1_FECollection(order, dim);
  FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec, dim, Ordering::byVDIM);
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


  ofstream mout("mass_mat.dat");
  m->PrintMatlab(mout);

  cout << "M.nnz = " << M.NumNonZeroElems() << endl;


  // Define a simple symmetric Gauss-Seidel preconditioner.
  GSSmoother prec(M);

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

    // Solve the SLAE
    const int print_iter = 0;
    PCG(M, prec, RHS, x_0, print_iter, 200, 1e-12, 0.0);

    // Compute and print the L^2 norm of the error
    if (time_step % tenth == 0)
      cout << "step " << time_step << " / " << n_time_steps
           << " ||solution||_{L^2} = " << x_0.Norml2() << endl;

    Vector nodal_values_x, nodal_values_y;
    x_0.GetNodalValues(nodal_values_x, 1);
    x_0.GetNodalValues(nodal_values_y, 2);

    if (visualization && (time_step % STEP_SNAP == 0))
      save_vts(time_step, "displ", sx, sy, nx, ny, nodal_values_x, nodal_values_y);

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

  delete b;
  delete m;
  delete s;
  delete fespace;
  delete fec;
}

