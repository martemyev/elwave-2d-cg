#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "config.hpp"
#include "mfem.hpp"

#include <string>
#include <vector>

static const char* DEFAULT_FILE_NAME = "no-file";

class ReceiversSet;



/**
 * Parameters descirbing the domain and the grid.
 */
class GridParameters
{
public:
  GridParameters();
  ~GridParameters() { }

  double sx, sy; ///< size of the computational domain
  int nx, ny; ///< number of cells in x- and y-directions

  void AddOptions(mfem::OptionsParser& args);
  void check_parameters() const;

private:
  GridParameters(const GridParameters&);
  GridParameters& operator=(const GridParameters&);
};



/**
 * Parameters describing the source.
 */
class SourceParameters
{
public:
  SourceParameters();
  ~SourceParameters() { }

  mfem::Vertex location;
  double frequency;
  double direction; ///< The direction of the point force source: 1 OX, 2 OY
  double scale;
  double Mxx, Mxy, Myy; ///< components of a moment tensor
  const char *type; ///< "pointforce", "momenttensor"
  const char *spatial_function; ///< "delta", "gauss"
  double gauss_support; ///< size of the support for the "gauss" spatial function
  bool plane_wave; ///< plane wave as a source at the depth of y-coordinate of
                   ///< the source location

  void AddOptions(mfem::OptionsParser& args);
  void check_parameters() const;

private:
  SourceParameters(const SourceParameters&);
  SourceParameters& operator=(const SourceParameters&);
};



/**
 * Parameters describing the media properties.
 */
class MediaPropertiesParameters
{
public:
  MediaPropertiesParameters();
  ~MediaPropertiesParameters();

  double rho, vp, vs; ///< homogeneous media properties

  const char *rhofile; ///< file names for heterogeneous media properties
  const char *vpfile;
  const char *vsfile;

  double *rho_array, *vp_array, *vs_array; ///< arrays of values describing
                                           ///< media properties

  double min_rho, max_rho, min_vp, max_vp, min_vs, max_vs;

  void AddOptions(mfem::OptionsParser& args);
  void check_parameters() const;
  void init(const GridParameters& grid);

private:
  MediaPropertiesParameters(const MediaPropertiesParameters&);
  MediaPropertiesParameters& operator=(const MediaPropertiesParameters&);
};



/**
 * Parameters describing the boundary conditions.
 */
class BoundaryConditionsParameters
{
public:
  BoundaryConditionsParameters();
  ~BoundaryConditionsParameters() { }

  const char* left; ///< left surface: absorbing (abs) or free
  const char* right; ///< right surface: absorbing (abs) or free
  const char* bottom; ///< bottom surface: absorbing (abs) or free
  const char* top; ///< top surface: absorbing (abs) or free
  double damp_layer; ///< thickness of a damping layer
  double damp_power; ///< power in damping coefficient functions

  void AddOptions(mfem::OptionsParser& args);
  void check_parameters() const;

private:
  BoundaryConditionsParameters(const BoundaryConditionsParameters&);
  BoundaryConditionsParameters& operator=(const BoundaryConditionsParameters&);
};



/**
 * Parameters of the problem to be solved.
 */
class Parameters
{
public:
  Parameters();
  ~Parameters();

  GridParameters grid;
  SourceParameters source;
  MediaPropertiesParameters media;
  BoundaryConditionsParameters bc;

  double T; ///< simulation time
  double dt; ///< time step
  int order; ///< finite element order

  int step_snap; ///< time step for outputting snapshots (every *th time step)
  const char* snapshot_format; ///< binary (bin) or VTK-type (vts)

  const char* method; ///< finite elements (fem) or spectral elements (sem)
  const char *extra_string; ///< added to output files for distinguishing the
                            ///< results

  const char *receivers_file; ///< file describing the sets of receivers
  std::vector<ReceiversSet*> sets_of_receivers;

  void init(int argc, char **argv);
  void check_parameters() const;

private:
  Parameters(const Parameters&); // no copies
  Parameters& operator=(const Parameters&); // no copies
};

#endif // PARAMETERS_HPP
