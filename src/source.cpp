#include "source.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;



Source::Source()
  : location(500.0, 500.0)
  , frequency(10.0)
  , angle(270.0)
  , ricker_scale(1.0)
  , point_force_scale(1.0)
  , gauss_support(1.0)
  , Mxx(1.0), Mxy(0.0), Myy(1.0) // explosive source
  , type(1)
{ }



void Source::AddOptions(OptionsParser& args)
{
  args.AddOption(&location(0), "-srcx", "--source-x", "x-coord of a source location");
  args.AddOption(&location(1), "-srcy", "--source-y", "y-coord of a source location");
  args.AddOption(&frequency, "-f0", "--frequency", "Central frequency of a source");
  args.AddOption(&angle, "-angle", "--angle", "Angle of a source force vector in degrees (0 is horizontal)");
  args.AddOption(&ricker_scale, "-rs", "--ricker-scale", "Factor for the Ricker wavelet");
  args.AddOption(&point_force_scale, "-pfs", "--point-force-scale", "Factor for the point force term of a source");
  args.AddOption(&gauss_support, "-gs", "--gauss-support", "Gauss support for Gaussian space functions of a source");
  args.AddOption(&Mxx, "-mxx", "--moment-tensor-xx", "xx-component of a moment tensor source");
  args.AddOption(&Mxy, "-mxy", "--moment-tensor-xy", "xy-component of a moment tensor source");
  args.AddOption(&Myy, "-myy", "--moment-tensor-yy", "yy-component of a moment tensor source");
  args.AddOption(&type, "-st", "--source-type", "Type of spatial source distribution (0 delta, 1 gauss)");
}



double Source::Ricker(double t) const
{
  const double a  = M_PI*frequency*(t-1./frequency);
  return ricker_scale * (1. - 2.*a*a)*exp(-a*a);
}



double Source::GaussFirstDerivative(double t) const
{
  const double a = M_PI*frequency*(t-1./frequency);
  return ricker_scale * (t-1./frequency)*exp(-a*a);
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

