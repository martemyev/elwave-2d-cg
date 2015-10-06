#include "parameters.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
Source::Source()
  : location(500.0, 500.0)
  , frequency(10.0)
  , angle(270.0)
  , ricker_scale(1.0)
  , point_force_scale(1.0)
  , gauss_support(1.0)
  , Mxx(1.0), Mxy(0.0), Myy(1.0)
  , type(1)
{ }

double Source::Ricker(double t) const
{
  const double a  = M_PI*frequency*(t-1./frequency);
  const double a2 = a*a;
  return ricker_scale * (1. - 2.*a2)*exp(-a2);
}

void Source::PointForce(const Vector& x, Vector& f) const
{
  if (type == 0) // Delta
    DeltaPointForce(x, f);
  else if (type == 1) // Gauss
    GaussPointForce(x, f);
  else
    mfem_error("Unknown source type");
}

void Source::MomentTensorSource(const Vector &x, Vector &f) const
{
  if (type == 0) // Delta
    DivDeltaMomentTensor(x, f);
  else if (type == 1) // Gauss
    DivGaussMomentTensor(x, f);
  else
    mfem_error("Unknown source type");
}

void Source::DeltaPointForce(const Vector& x, Vector& f) const
{
  const double tol = 1e-2;
  const double loc[] = { location(0), location(1) };
  if (x.DistanceTo(loc) < tol)
  {
    f(0) = point_force_scale * cos(to_radians(angle));
    f(1) = point_force_scale * sin(to_radians(angle));
  }
  else
    f = 0.0;
}

void Source::GaussPointForce(const Vector& x, Vector& f) const
{
  const double xdiff  = x(0)-location(0);
  const double ydiff  = x(1)-location(1);
  const double xdiff2 = xdiff*xdiff;
  const double ydiff2 = ydiff*ydiff;
  const double h2 = gauss_support*gauss_support;
  const double G = exp(-(xdiff2 + ydiff2) / h2);
  f(0) = point_force_scale * G*cos(to_radians(angle));
  f(1) = point_force_scale * G*sin(to_radians(angle));
}

void Source::DivDeltaMomentTensor(const Vector& x, Vector& f) const
{
  mfem_error("NOT implemented");
}

void Source::DivGaussMomentTensor(const Vector& x, Vector& f) const
{
  const double xdiff  = x(0)-location(0);
  const double ydiff  = x(1)-location(1);
  const double xdiff2 = xdiff*xdiff;
  const double ydiff2 = ydiff*ydiff;
  const double h2 = gauss_support*gauss_support;
  const double exp_val = exp(-(xdiff2 + ydiff2) / h2);
  const double Gx = -2.*xdiff/h2 * exp_val;
  const double Gy = -2.*ydiff/h2 * exp_val;

  f(0) = Mxx*Gx + Mxy*Gy;
  f(1) = Mxy*Gx + Myy*Gy;
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
Parameters::Parameters()
  : sx(1000.0)
  , sy(1000.0)
  , nx(10)
  , ny(10)
  , T(1.0)
  , dt(1e-3)
  , order(1)
  , rho(1000)
  , vp(3000)
  , vs(1700)
  , rhofile(DEFAULT_FILE_NAME)
  , vpfile(DEFAULT_FILE_NAME)
  , vsfile(DEFAULT_FILE_NAME)
  , rho_array(nullptr)
  , vp_array(nullptr)
  , vs_array(nullptr)
  , damp_layer(100.0)
  , topsurf(1)
  , source()
  , step_snap(1000)
{ }

Parameters::~Parameters()
{
  delete[] vs_array;
  delete[] vp_array;
  delete[] rho_array;
}

void Parameters::init(int argc, char **argv)
{
  OptionsParser args(argc, argv);
  args.AddOption(&sx, "-sx", "--sizex", "Size of domain in x-direction, m");
  args.AddOption(&sy, "-sy", "--sizey", "Size of domain in y-direction, m");
  args.AddOption(&nx, "-nx", "--numberx", "Number of elements in x-direction");
  args.AddOption(&ny, "-ny", "--numbery", "Number of elements in y-direction");
  args.AddOption(&T, "-T", "--time-end", "Simulation time, s");
  args.AddOption(&dt, "-dt", "--time-step", "Time step, s");
  args.AddOption(&order, "-o", "--order", "Finite element order (polynomial degree)");

  args.AddOption(&rho, "-rho", "--rho", "Density of homogeneous model, kg/m^3");
  args.AddOption(&vp, "-vp", "--vp", "P-wave velocity of homogeneous model, m/s");
  args.AddOption(&vs, "-vs", "--vs", "S-wave velocity of homogeneous model, m/s");
  args.AddOption(&rhofile, "-rhofile", "--rhofile", "Density file, in kg/m^3");
  args.AddOption(&vpfile, "-vpfile", "--vpfile", "P-wave velocity file, in m/s");
  args.AddOption(&vsfile, "-vsfile", "--vsfile", "S-wave velocity file, in m/s");

  args.AddOption(&damp_layer, "-damp", "--damp-layer", "Thickness of damping layer, m");
  args.AddOption(&topsurf, "-top", "--top-surface", "Top surface: 0 absorbing, 1 free");

  args.AddOption(&source.location(0), "-srcx", "--source-x", "x-coord of a source location");
  args.AddOption(&source.location(1), "-srcy", "--source-y", "y-coord of a source location");
  args.AddOption(&source.frequency, "-f0", "--frequency", "Central frequency of a source");
  args.AddOption(&source.angle, "-angle", "--angle", "Angle of a source force vector in degrees (0 is horizontal)");
  args.AddOption(&source.ricker_scale, "-rs", "--ricker-scale", "Factor for the Ricker wavelet");
  args.AddOption(&source.point_force_scale, "-pfs", "--point-force-scale", "Factor for the point force term of a source");
  args.AddOption(&source.gauss_support, "-gs", "--gauss-support", "Gauss support for Gaussian space functions of a source");
  args.AddOption(&source.Mxx, "-mxx", "--moment-tensor-xx", "xx-component of a moment tensor source");
  args.AddOption(&source.Mxy, "-mxy", "--moment-tensor-xy", "xy-component of a moment tensor source");
  args.AddOption(&source.Myy, "-myy", "--moment-tensor-yy", "yy-component of a moment tensor source");
  args.AddOption(&source.type, "-st", "--source-type", "Type of spatial source distribution (0 delta, 1 gauss)");

  args.AddOption(&step_snap, "-step-snap", "--step-snap", "Time step for outputting snapshots");

  args.Parse();
  if (!args.Good())
  {
    args.PrintUsage(cout);
    throw 1;
  }
  args.PrintOptions(cout);

  const int n_elements = nx*ny;

  rho_array = new double[n_elements];
  vp_array = new double[n_elements];
  vs_array = new double[n_elements];

  double min_rho, max_rho, min_vp, max_vp, min_vs, max_vs;

  if (!strcmp(rhofile, DEFAULT_FILE_NAME))
  {
    for (int i = 0; i < n_elements; ++i) rho_array[i] = rho;
    min_rho = max_rho = rho;
  }
  else
  {
    read_binary(rhofile, n_elements, rho_array);
    get_minmax(rho_array, n_elements, min_rho, max_rho);
  }

  if (!strcmp(vpfile, DEFAULT_FILE_NAME))
  {
    for (int i = 0; i < n_elements; ++i) vp_array[i] = vp;
    min_vp = max_vp = vp;
  }
  else
  {
    read_binary(vpfile, n_elements, vp_array);
    get_minmax(vp_array, n_elements, min_vp, max_vp);
  }

  if (!strcmp(vsfile, DEFAULT_FILE_NAME))
  {
    for (int i = 0; i < n_elements; ++i) vs_array[i] = vs;
    min_vs = max_vs = vs;
  }
  else
  {
    read_binary(vsfile, n_elements, vs_array);
    get_minmax(vs_array, n_elements, min_vs, max_vs);
  }

  const double min_wavelength = min(min_vp, min_vs) / (2.0*source.frequency);
  cout << "min wavelength = " << min_wavelength << endl;

  if (damp_layer < 2.5*min_wavelength)
    mfem_warning("damping layer for absorbing bc should be about 3*wavelength");
}