#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "config.hpp"
#include "mfem.hpp"

#include <string>
#include <vector>

static const char* DEFAULT_FILE_NAME = "no-file";

class ReceiversSet;


class Source
{
public:
  Source();

  mfem::Vertex location;
  double frequency;
  double angle;
  double ricker_scale;
  double point_force_scale;
  double gauss_support;
  double Mxx, Mxy, Myy; // components of a moment tensor
  int type; // 0 - Delta function, 1 - Gaussian function

  void AddOptions(mfem::OptionsParser& args);

  double Ricker(double t) const;
  void PointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void MomentTensorSource(const mfem::Vector& x, mfem::Vector& f) const;

private:
  void DeltaPointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void GaussPointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void DivDeltaMomentTensor(const mfem::Vector& x, mfem::Vector& f) const;
  void DivGaussMomentTensor(const mfem::Vector& x, mfem::Vector& f) const;
};



class Parameters
{
public:
  Parameters();
  ~Parameters();

  void init(int argc, char **argv);

  double sx, sy; // size of the computational domain
  int nx, ny; // number of cells in x- and y-directions
  double T; // simulation time
  double dt; // time step
  int order; // finite element order

  double rho, vp, vs; // homogeneous media properties

  const char *rhofile; // heterogeneous media properties
  const char *vpfile;
  const char *vsfile;

  // arrays of values describing media properties
  double *rho_array, *vp_array, *vs_array;

  double damp_layer; // thickness of a damping layer
  int topsurf; // top surface: 0 absorbing, 1 free

  Source source; // source of the wave

  int step_snap; // time step for outputting snapshots (every *th time step)

  int method; // 0 - FEM, 1 - SEM

  const char *extra_string;

  const char *receivers_file;
  std::vector<ReceiversSet*> sets_of_receivers;


private:
  Parameters(const Parameters&); // no copies
  Parameters& operator=(const Parameters&); // no copies
};

#endif // PARAMETERS_HPP
