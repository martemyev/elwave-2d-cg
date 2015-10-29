#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "config.hpp"
#include "mfem.hpp"
#include "source.hpp"

#include <string>
#include <vector>

static const char* DEFAULT_FILE_NAME = "no-file";

class ReceiversSet;



class Parameters
{
public:
  Parameters();
  ~Parameters();

  void init(int argc, char **argv);

  double sx, sy; ///< size of the computational domain
  int nx, ny; ///< number of cells in x- and y-directions
  double T; ///< simulation time
  double dt; ///< time step
  int order; ///< finite element order

  double rho, vp, vs; ///< homogeneous media properties

  const char *rhofile; ///< heterogeneous media properties
  const char *vpfile;
  const char *vsfile;

  double *rho_array, *vp_array, *vs_array; ///< arrays of values describing
                                           ///< media properties

  double damp_layer; ///< thickness of a damping layer
  double damp_power; ///< power in damping coefficient functions
  const char* topsurf; ///< top surface: absorbing (abs) or free

  Source source; ///< source of the wave

  int step_snap; ///< time step for outputting snapshots (every *th time step)
  const char* snapshot_format; ///< binary (bin) or VTK-type (vts)

  const char* method; ///< finite elements (fem) or spectral elements (sem)

  const char *extra_string; ///< added to output files for distinguishing the
                            ///< results

  const char *receivers_file; ///< file describing the sets of receivers
  std::vector<ReceiversSet*> sets_of_receivers;


private:
  Parameters(const Parameters&); // no copies
  Parameters& operator=(const Parameters&); // no copies
};

#endif // PARAMETERS_HPP
