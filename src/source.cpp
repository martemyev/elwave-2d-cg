#include "source.hpp"
#include "utilities.hpp"
#include "parameters.hpp"

using namespace std;
using namespace mfem;



Source::Source()
  : location(500.0, 500.0)
  , frequency(10.0)
  , direction(2) // OY
  , scale(1e+6)
  , Mxx(1.0), Mxy(0.0), Myy(1.0) // explosive source
  , type_string("pointforce")
  , spatial_function("gauss")
  , gauss_support(10.0)
  , plane_wave(false)
  , type(POINT_FORCE)
{ }



void Source::AddOptions(OptionsParser& args)
{
  args.AddOption(&location(0), "-srcx", "--source-x", "x-coord of a source location");
  args.AddOption(&location(1), "-srcy", "--source-y", "y-coord of a source location");
  args.AddOption(&frequency, "-f0", "--source-frequency", "Central frequency of a source");
  args.AddOption(&direction, "-dir", "--source-direction", "Direction of the point force source (1 OX, 2 OY)");
  args.AddOption(&scale, "-scale", "--source-scale", "Scaling factor for the source");
  args.AddOption(&Mxx, "-mxx", "--moment-tensor-xx", "xx-component of a moment tensor source");
  args.AddOption(&Mxy, "-mxy", "--moment-tensor-xy", "xy-component of a moment tensor source");
  args.AddOption(&Myy, "-myy", "--moment-tensor-yy", "yy-component of a moment tensor source");
  args.AddOption(&type_string, "-type", "--source-type", "Type of the source (pointforce, momenttensor)");
  args.AddOption(&spatial_function, "-spatial", "--source-spatial", "Spatial function of the source (delta, gauss)");
  args.AddOption(&gauss_support, "-gs", "--gauss-support", "Gauss support for 'gauss' spatial function of the source");
  args.AddOption(&plane_wave, "-planewave", "--plane-wave", "-noplanewave", "--no-plane-wave", "Plane wave as a source");
}



void Source::check_and_update_parameters()
{
  MFEM_VERIFY(frequency > 0, "Frequency (" + d2s(frequency) + ") must be >0");
  MFEM_VERIFY(direction == 1 || direction == 2, "Unsupported direction of the "
              "source: " + d2s(direction));

  if (!strcmp(type_string, "pointforce"))
    type = POINT_FORCE;
  else if (!strcmp(type_string, "momenttensor"))
    type = MOMENT_TENSOR;
  else
    MFEM_ABORT("Unknown source type: " + string(type_string));

  MFEM_VERIFY(!strcmp(spatial_function, "delta") ||
              !strcmp(spatial_function, "gauss"), "Unknown spatial function of "
              "the source: " + string(spatial_function));
}



double Source::Ricker(double t) const
{
  const double a  = M_PI*frequency*(t-1./frequency);
  return scale*(1.-2.*a*a)*exp(-a*a);
}



double Source::GaussFirstDerivative(double t) const
{
  const double a = M_PI*frequency*(t-1./frequency);
  return scale*(t-1./frequency)*exp(-a*a);
}



void Source::PointForce(const Vector& source_location, const Vector& x,
                        Vector& f) const
{
  if (!strcmp(spatial_function, "delta"))
    DeltaPointForce(source_location, x, f);
  else if (!strcmp(spatial_function, "gauss"))
    GaussPointForce(source_location, x, f);
  else
    MFEM_ABORT("Unknown source spatial function: " + string(spatial_function));
}



void Source::MomentTensorSource(const Vector& source_location, const Vector &x,
                                Vector &f) const
{
  if (!strcmp(spatial_function, "delta"))
    DivDeltaMomentTensor(source_location, x, f);
  else if (!strcmp(spatial_function, "gauss"))
    DivGaussMomentTensor(source_location, x, f);
  else
    MFEM_ABORT("Unknown source spatial function: " + string(spatial_function));
}



void Source::PlaneWaveSource(const Parameters& param, const Vector &x,
                             Vector &f) const
{
  const double px = x(0); // coordinates of the point
  const double py = x(1);

  const double X0 = 0.0;
  const double X1 = param.sx;
  const double layer = param.damp_layer;
  const double tol = FLOAT_NUMBERS_EQUALITY_TOLERANCE;

  // if the point 'x' is on the plane of the plane wave, and it's not in the
  // absorbing layer, we have a source located at the exact same point
  if (fabs(py - location(1)) < tol) // &&
      //px-layer+tol > X0 && px+layer-tol < X1)
  {
    if (type == POINT_FORCE)
      PointForce(x, x, f);
    else if (type == MOMENT_TENSOR)
      MomentTensorSource(x, x, f);
    else MFEM_ABORT("Unknown source type: " + string(type_string));
  }
  else
    f = 0.0;
}



void Source::DeltaPointForce(const Vector& source_location, const Vector& x,
                             Vector& f) const
{
  const double tol = FLOAT_NUMBERS_EQUALITY_TOLERANCE;
  const double loc[] = { source_location(0), source_location(1) };
  double value = 0.0;
  if (x.DistanceTo(loc) < tol)
    value = 1.0;

  f = 0.0;
  f(direction-1) = value;
}



void Source::GaussPointForce(const Vector& source_location, const Vector& x,
                             Vector& f) const
{
  const double xdiff  = x(0) - source_location(0);
  const double ydiff  = x(1) - source_location(1);
  const double xdiff2 = xdiff*xdiff;
  const double ydiff2 = ydiff*ydiff;
  const double h2 = gauss_support*gauss_support;
  const double G = exp(-(xdiff2 + ydiff2) / h2);
  f = 0.0;
  f(direction-1) = G;
}



void Source::DivDeltaMomentTensor(const Vector&, const Vector&, Vector&) const
{
  mfem_error("NOT implemented");
}



void Source::DivGaussMomentTensor(const Vector& source_location,
                                  const Vector& x, Vector& f) const
{
  const double xdiff  = x(0) - source_location(0);
  const double ydiff  = x(1) - source_location(1);
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
VectorPointForce::VectorPointForce(int dim, const Source& s)
  : VectorCoefficient(dim)
  , source(s)
{
  location.SetSize(vdim);
  for (int i = 0; i < vdim; ++i)
    location(i) = source.location(i);
}



void VectorPointForce::Eval(Vector &V, ElementTransformation &T,
                            const IntegrationPoint &ip)
{
  Vector transip;
  T.Transform(ip, transip);
  V.SetSize(vdim);
  Vector location(vdim);
  for (int i = 0; i < vdim; ++i)
    location(i) = source.location(i);
  source.PointForce(location, transip, V);
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
MomentTensorSource::MomentTensorSource(int dim, const Source& s)
  : VectorCoefficient(dim)
  , source(s)
{
  location.SetSize(vdim);
  for (int i = 0; i < vdim; ++i)
    location(i) = source.location(i);
}



void MomentTensorSource::Eval(Vector &V, ElementTransformation &T,
                              const IntegrationPoint &ip)
{
  Vector transip;
  T.Transform(ip, transip);
  V.SetSize(vdim);
  source.MomentTensorSource(location, transip, V);
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
PlaneWaveSource::PlaneWaveSource(int dim, const Parameters& p)
  : VectorCoefficient(dim)
  , param(p)
{ }



void PlaneWaveSource::Eval(Vector &V, ElementTransformation &T,
                           const IntegrationPoint &ip)
{
  Vector transip;
  T.Transform(ip, transip);
  V.SetSize(vdim);
  param.source.PlaneWaveSource(param, transip, V);
}
