#include "parameters.hpp"
#include "receivers.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cfloat>

using namespace std;
using namespace mfem;



//------------------------------------------------------------------------------
//
// Grid parameters
//
//------------------------------------------------------------------------------
GridParameters::GridParameters()
  : sx(1000.0)
  , sy(1000.0)
  , nx(10)
  , ny(10)
{ }

void GridParameters::AddOptions(OptionsParser& args)
{
  args.AddOption(&sx, "-sx", "--sizex", "Size of domain in x-direction, m");
  args.AddOption(&sy, "-sy", "--sizey", "Size of domain in y-direction, m");
  args.AddOption(&nx, "-nx", "--numberx", "Number of elements in x-direction");
  args.AddOption(&ny, "-ny", "--numbery", "Number of elements in y-direction");
}

void GridParameters::check_parameters() const
{
  MFEM_VERIFY(sx > 0 && sy > 0, "Size of the domain (sx=" + d2s(sx) + " m, sy="+
              d2s(sy) + " m) must be >0");
  MFEM_VERIFY(nx > 0 && ny > 0, "Number of cells (nx=" + d2s(nx) + ", ny=" +
              d2s(ny) + ") must be >0");
}



//------------------------------------------------------------------------------
//
// Source parameters
//
//------------------------------------------------------------------------------
SourceParameters::SourceParameters()
  : location(500.0, 500.0)
  , frequency(10.0)
  , direction(2) // OY
  , scale(1e+6)
  , Mxx(1.0), Mxy(0.0), Myy(1.0) // explosive source
  , type("pointforce")
  , spatial_function("gauss")
  , time_function("auto")
  , gauss_support(10.0)
  , plane_wave(false)
{ }

void SourceParameters::AddOptions(OptionsParser& args)
{
  args.AddOption(&location(0), "-srcx", "--source-x", "x-coord of a source location");
  args.AddOption(&location(1), "-srcy", "--source-y", "y-coord of a source location");
  args.AddOption(&frequency, "-f0", "--source-frequency", "Central frequency of a source");
  args.AddOption(&direction, "-dir", "--source-direction", "Direction of the point force source (1 OX, 2 OY)");
  args.AddOption(&scale, "-scale", "--source-scale", "Scaling factor for the source");
  args.AddOption(&Mxx, "-mxx", "--moment-tensor-xx", "xx-component of a moment tensor source");
  args.AddOption(&Mxy, "-mxy", "--moment-tensor-xy", "xy-component of a moment tensor source");
  args.AddOption(&Myy, "-myy", "--moment-tensor-yy", "yy-component of a moment tensor source");
  args.AddOption(&type, "-type", "--source-type", "Type of the source (pointforce, momenttensor)");
  args.AddOption(&spatial_function, "-spatial", "--source-spatial", "Spatial function of the source (delta, gauss)");
  args.AddOption(&time_function, "-temporary", "--source-temporary", "Temporary function of the source (auto, ricker, firstgauss, gauss)");
  args.AddOption(&gauss_support, "-gs", "--gauss-support", "Gauss support for 'gauss' spatial function of the source");
  args.AddOption(&plane_wave, "-planewave", "--plane-wave", "-noplanewave", "--no-plane-wave", "Plane wave as a source");
}

void SourceParameters::check_parameters() const
{
  MFEM_VERIFY(frequency > 0, "Frequency (" + d2s(frequency) + ") must be >0");
  MFEM_VERIFY(direction == 1 || direction == 2, "Unsupported direction of the "
              "source: " + d2s(direction));
  MFEM_VERIFY(!strcmp(type, "pointforce") || !strcmp(type, "momenttensor"),
              "Unknown source type: " + string(type));
  MFEM_VERIFY(!strcmp(spatial_function, "delta") ||
              !strcmp(spatial_function, "gauss"), "Unknown spatial function of "
              "the source: " + string(spatial_function));
  if (!strcmp(spatial_function, "gauss"))
    MFEM_VERIFY(gauss_support > 0, "Gauss support (" + d2s(gauss_support) +
                ") must be >0");

  MFEM_VERIFY(!strcmp(time_function, "auto") ||
              !strcmp(time_function, "ricker") ||
              !strcmp(time_function, "firstgauss") ||
              !strcmp(time_function, "gauss"), "Unknown temporary function "
              "of the source: " + string(time_function));
}



//------------------------------------------------------------------------------
//
// Media properties parameters
//
//------------------------------------------------------------------------------
MediaPropertiesParameters::MediaPropertiesParameters()
  : rho(2500.0)
  , vp(3500.0)
  , vs(2000.0)
  , rhofile(DEFAULT_FILE_NAME)
  , vpfile(DEFAULT_FILE_NAME)
  , vsfile(DEFAULT_FILE_NAME)
  , rho_array(nullptr)
  , vp_array(nullptr)
  , vs_array(nullptr)
  , min_rho(DBL_MAX), max_rho(DBL_MIN)
  , min_vp (DBL_MAX), max_vp (DBL_MIN)
  , min_vs (DBL_MAX), max_vs (DBL_MIN)
{ }

MediaPropertiesParameters::~MediaPropertiesParameters()
{
  delete[] rho_array;
  delete[] vp_array;
  delete[] vs_array;
}

void MediaPropertiesParameters::AddOptions(OptionsParser& args)
{
  args.AddOption(&rho, "-rho", "--rho", "Density of homogeneous model, kg/m^3");
  args.AddOption(&vp,  "-vp",  "--vp",  "P-wave velocity of homogeneous model, m/s");
  args.AddOption(&vs,  "-vs",  "--vs",  "S-wave velocity of homogeneous model, m/s");
  args.AddOption(&rhofile, "-rhofile", "--rhofile", "Density file, in kg/m^3");
  args.AddOption(&vpfile,  "-vpfile",  "--vpfile",  "P-wave velocity file, in m/s");
  args.AddOption(&vsfile,  "-vsfile",  "--vsfile",  "S-wave velocity file, in m/s");
}

void MediaPropertiesParameters::check_parameters() const
{
  // no checks here
}

void MediaPropertiesParameters::init(const GridParameters& grid)
{
  const int n_elements = grid.nx*grid.ny;

  rho_array = new double[n_elements];
  vp_array = new double[n_elements];
  vs_array = new double[n_elements];

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
}



//------------------------------------------------------------------------------
//
// Boundary conditions parameters
//
//------------------------------------------------------------------------------
BoundaryConditionsParameters::BoundaryConditionsParameters()
  : left("abs")
  , right("abs")
  , bottom("abs")
  , top("free")
  , damp_layer(100.0)
  , damp_power(3.0)
{ }

void BoundaryConditionsParameters::AddOptions(OptionsParser& args)
{
  // Left, right and bottom surfaces are usually absorbing, so we don't need to
  // set up program options for them, but this can be changed if desired.
//  args.AddOption(&left, "-left", "--left-surface", "Left surface: abs or free");
//  args.AddOption(&right, "-right", "--right-surface", "Right surface: abs or free");
//  args.AddOption(&bottom, "-bottom", "--bottom-surface", "Bottom surface: abs or free");

  args.AddOption(&top, "-top", "--top-surface", "Top surface: abs or free");
  args.AddOption(&damp_layer, "-dlayer", "--damp-layer", "Thickness of damping layer, m");
  args.AddOption(&damp_power, "-dpower", "--damp-power", "Power in damping coefficient functions");
}

void BoundaryConditionsParameters::check_parameters() const
{
  MFEM_VERIFY(!strcmp(left, "abs") || !strcmp(left, "free"), "Unknown boundary "
              "condition on the left surface: " + string(left));
  MFEM_VERIFY(!strcmp(right, "abs") || !strcmp(right, "free"), "Unknown boundary "
              "condition on the right surface: " + string(right));
  MFEM_VERIFY(!strcmp(bottom, "abs") || !strcmp(bottom, "free"), "Unknown boundary "
              "condition on the bottom surface: " + string(bottom));
  MFEM_VERIFY(!strcmp(top, "abs") || !strcmp(top, "free"), "Unknown boundary "
              "condition on the top surface: " + string(top));
  if (!strcmp(left, "abs") || !strcmp(right, "abs") || !strcmp(bottom, "abs") ||
      !strcmp(top, "abs"))
    MFEM_VERIFY(damp_layer > 0, "Damping layer (" + d2s(damp_layer) +
                ") must be >0");
}



//------------------------------------------------------------------------------
//
// Reservoir parameters
//
//------------------------------------------------------------------------------
ReservoirParameters::ReservoirParameters()
  : exists(false)                     // no reservoir by default
  , nx(0), ny(0)                      // if not set up - no reservoir
  , x0(0.0), x1(0.0)
  , y0(0.0), y1(0.0)
  , rho(2500.0)
  , vp(3500.0)
  , vs(2000.0)
  , rhofile(DEFAULT_FILE_NAME)
  , vpfile(DEFAULT_FILE_NAME)
  , vsfile(DEFAULT_FILE_NAME)
  , extend_to_damp(false)
  , rho_array(nullptr)
  , vp_array(nullptr)
  , vs_array(nullptr)
  , min_rho(DBL_MAX), max_rho(DBL_MIN)
  , min_vp (DBL_MAX), max_vp (DBL_MIN)
  , min_vs (DBL_MAX), max_vs (DBL_MIN)
{ }

ReservoirParameters::~ReservoirParameters()
{
  delete[] rho_array;
  delete[] vp_array;
  delete[] vs_array;
}

void ReservoirParameters::AddOptions(OptionsParser &args)
{
  args.AddOption(&nx, "-rsr-nx", "--reservoir-nx", "Number of cells in x-direction in reservoir layer");
  args.AddOption(&ny, "-rsr-ny", "--reservoir-ny", "Number of cells in y-direction in reservoir layer");
  args.AddOption(&x0, "-rsr-x0", "--reservoir-x0", "Left boundary of the reservoir layer, m");
  args.AddOption(&x1, "-rsr-x1", "--reservoir-x1", "Right boundary of the reservoir layer, m");
  args.AddOption(&y0, "-rsr-y0", "--reservoir-y0", "Bottom boundary of the reservoir layer, m");
  args.AddOption(&y1, "-rsr-y1", "--reservoir-y1", "top boundary of the reservoir layer, m");
  args.AddOption(&rho, "-rsr-rho", "--reservoir-rho", "Density of homogeneous reservoir, kg/m^3");
  args.AddOption(&vp,  "-rsr-vp",  "--reservoir-vp",  "P-wave velocity of homogeneous reservoir, m/s");
  args.AddOption(&vs,  "-rsr-vs",  "--reservoir-vs",  "S-wave velocity of homogeneous reservoir, m/s");
  args.AddOption(&rhofile, "-rsr-rhofile", "--reservoir-rhofile", "Density file for reservoir, in kg/m^3");
  args.AddOption(&vpfile,  "-rsr-vpfile",  "--reservoir-vpfile",  "P-wave velocity file for reservoir, in m/s");
  args.AddOption(&vsfile,  "-rsr-vsfile",  "--reservoir-vsfile",  "S-wave velocity file for reservoir, in m/s");
  args.AddOption(&extend_to_damp,  "-rsr-ext-damp", "--reservoir-extend-damp",
                 "-no-rsr-ext-damp", "--no-reservoir-extend-damp",
                 "Extend reservoir layer so it covers damping layer");
}

void ReservoirParameters::check_parameters(const GridParameters& grid) const
{
  if (nx*ny > 0)
  {
    MFEM_VERIFY(x0 > 0, "Left boundary of the reservoir (" + d2s(x0) +
                ") should be >0");
    MFEM_VERIFY(x1 < grid.sx, "Right boundary of the reservoir (" + d2s(x1) +
                ") should be <sx (size domain in x-direction: " +
                d2s(grid.sx) + ")");
    MFEM_VERIFY(x1 > x0, "x1 should be >x0 for reservoir");
    MFEM_VERIFY(y0 > 0, "Bottom boundary of the reservoir (" + d2s(y0) +
                ") should be >0");
    MFEM_VERIFY(y1 < grid.sy, "Top boundary of the reservoir (" + d2s(y1) +
                ") should be <sy (size domain in y-direction: " +
                d2s(grid.sy) + ")");
    MFEM_VERIFY(y1 > y0, "y1 should be >y0 for reservoir");
  }
}

void ReservoirParameters::init()
{
  if (nx*ny > 0) exists = true;
  else return; // if the reservoir doesn't exist - there is nothing to initialize

  const int n_elements = nx*ny;

  rho_array = new double[n_elements];
  vp_array  = new double[n_elements];
  vs_array  = new double[n_elements];

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
}

bool ReservoirParameters::contains(const Vector &point) const
{
  if (!exists) return false;
  const double tol = FIND_CELL_TOLERANCE;
  const double y = point(1);
  if (y+tol < y0 || y-tol > y1) return false; // not in the layer
  if (!extend_to_damp) // if we restrict the reservoir domain (it's not a layer)
  {                    // then we need to check x-coordinate too
    const double x = point(0);
    if (x+tol < x0 || x-tol > x1) return false;
  }
  // otherwise the point is in the reservoir
  return true;
}

int ReservoirParameters::find_cell(const Vector &point) const
{
  MFEM_ASSERT(exists, "There is no reservoir");
  const double tol = FIND_CELL_TOLERANCE;

  const double y = point(1);
  int celly = (y-y0)*ny/(y1-y0);
  if (celly) --celly;

  const double x = point(0);
  if (extend_to_damp)
  {
    if (x+tol < x0) return celly*nx;
    if (x-tol > x1) return celly*nx+nx-1;
  }

  // if the point is in the reservoir
  int cellx = (x-x0)*nx/(x1-x0);
  if (cellx) --cellx;

  return celly*nx+cellx;
}

string ReservoirParameters::info() const
{
  string str = "\nReservoir properties:\nnx = " + d2s(nx) + " ny = " + d2s(ny)+
               "\nx in [" + d2s(x0) + ", " + d2s(x1) + "]\ny in [" + d2s(y0) +
               ", " + d2s(y1) + "]\nrho: min = " + d2s(min_rho) + ", max = " +
               d2s(max_rho) + "\nvp: min = " + d2s(min_vp) + ", max = " +
               d2s(max_vp) + "\nvs: min = " + d2s(min_vs) + ", max = " +
               d2s(max_vs) + "\n";
  return str;
}



//------------------------------------------------------------------------------
//
// All parameters of the problem to be solved
//
//------------------------------------------------------------------------------
Parameters::Parameters()
  : grid()
  , source()
  , media()
  , bc()
  , T(1.0)
  , dt(1e-3)
  , order(1)
  , step_snap(1000)
  , snapshot_format("vts")
  , method("sem")
  , extra_string("")
  , receivers_file(DEFAULT_FILE_NAME)
{ }

Parameters::~Parameters()
{
  for (size_t i = 0; i < sets_of_receivers.size(); ++i)
    delete sets_of_receivers[i];
}

void Parameters::init(int argc, char **argv)
{
  OptionsParser args(argc, argv);

  grid.AddOptions(args);
  source.AddOptions(args);
  media.AddOptions(args);
  bc.AddOptions(args);
  reservoir.AddOptions(args);

  args.AddOption(&T, "-T", "--time-end", "Simulation time, s");
  args.AddOption(&dt, "-dt", "--time-step", "Time step, s");
  args.AddOption(&order, "-o", "--order", "Finite element order (polynomial degree)");
  args.AddOption(&step_snap, "-step-snap", "--step-snapshot", "Time step for outputting snapshots");
  args.AddOption(&snapshot_format, "-snap-format", "--snapshot-format", "Format of snapshots (bin, vts)");
  args.AddOption(&method, "-method", "--method", "Finite elements (fem) or spectral elements (sem)");
  args.AddOption(&extra_string, "-extra", "--extra", "Extra string for naming output files");
  args.AddOption(&receivers_file, "-rec-file", "--receivers-file", "File with information about receivers");

  args.Parse();
  if (!args.Good())
  {
    args.PrintUsage(cout);
    throw 1;
  }
  args.PrintOptions(cout);

  check_parameters();

  media.init(grid);
  reservoir.init();

  cout << reservoir.info() << endl;

  vector<double> min_velocities;
  min_velocities.push_back(media.min_vp);
  min_velocities.push_back(media.min_vs);
  if (reservoir.exists)
  {
    min_velocities.push_back(reservoir.min_vp);
    min_velocities.push_back(reservoir.min_vs);
  }
  const double v_min = *min_element(min_velocities.begin(), min_velocities.end());
  const double min_wavelength = v_min / (2.0*source.frequency);
  cout << "min wavelength = " << min_wavelength << endl;

  if (bc.damp_layer < 2.5*min_wavelength)
    mfem_warning("damping layer for absorbing bc should be about 3*wavelength");

  ifstream in(receivers_file);
  MFEM_VERIFY(in, "The file '" + string(receivers_file) + "' can't be opened");
  string line; // we read the file line-by-line
  string type; // type of the set of receivers
  while (getline(in, line))
  {
    // ignore empty lines and lines starting from '#'
    if (line.empty() || line[0] == '#') continue;
    // every meaningfull line should start with the type of the receivers set
    istringstream iss(line);
    iss >> type;
    ReceiversSet *rec_set = nullptr;
    if (type == "Line")
      rec_set = new ReceiversLine();
    else MFEM_ABORT("Unknown type of receivers set: " + type);

    rec_set->init(in); // read the parameters
    rec_set->distribute_receivers();
    rec_set->find_cells_containing_receivers(grid.nx, grid.ny, grid.sx, grid.sy);
    sets_of_receivers.push_back(rec_set); // put this set in the vector
  }
}

void Parameters::check_parameters() const
{
  grid.check_parameters();
  source.check_parameters();
  media.check_parameters();
  bc.check_parameters();
  reservoir.check_parameters(grid);

  MFEM_VERIFY(T > 0, "Time (" + d2s(T) + ") must be >0");
  MFEM_VERIFY(dt < T, "dt (" + d2s(dt) + ") must be < T (" + d2s(T) + ")");
  MFEM_VERIFY(order > 0, "FEM order (" + d2s(order) + ") must be >0");
  MFEM_VERIFY(step_snap > 0, "step_snap (" + d2s(step_snap) + ") must be >0");
  MFEM_VERIFY(!strcmp(method, "fem") || !strcmp(method, "FEM") ||
              !strcmp(method, "sem") || !strcmp(method, "SEM"), "Method (" +
              string(method) + ") must be either fem or sem");
}
