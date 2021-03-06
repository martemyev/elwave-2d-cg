#ifndef ELWAVE2D_HPP
#define ELWAVE2D_HPP

#include "config.hpp"
#include "mfem.hpp"

#include <vector>

class Parameters;



/**
 * 2D elastic wave run by finite elements or spectral elements (both are
 * continuous Galerkin approaches).
 */
class ElasticWave2D
{
public:
  ElasticWave2D(const Parameters& p) : param(p) { }
  ~ElasticWave2D() { }

  void run();

private:
  const Parameters& param;

  /**
   * Finite Element Method (FEM) (non-diagonal mass matrix) with Absorbing
   * Layers by Increasing Damping (ALID) for implementation of absorbing
   * boundary condition.
   */
  void run_FEM_ALID();

  /**
   * Spectral Element Method (SEM) (diagonal mass matrix) with Stiffness
   * Reduction Method (SRM) for implementation of absorbing boundary condition.
   */
  void run_SEM_SRM();
};



/**
 * Cell-wise constant coefficient.
 */
class CWConstCoefficient : public mfem::Coefficient
{
public:
  CWConstCoefficient(double *array, bool own = 1)
    : val_array(array), own_array(own)
  { }

  virtual ~CWConstCoefficient() { if (own_array) delete[] val_array; }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    return val_array[T.ElementNo];
  }

protected:
  double *val_array;
  bool own_array;
};



/**
 * A coefficient obtained with multiplication of a cell-wise constant
 * coefficient and a function.
 */
class CWFunctionCoefficient : public CWConstCoefficient
{
public:
  CWFunctionCoefficient(double(*F)(const mfem::Vector&, const Parameters&),
                        const Parameters& _param,
                        double *array, bool own = 1)
    : CWConstCoefficient(array, own)
    , Function(F)
    , param(_param)
  { }
  virtual ~CWFunctionCoefficient() { }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    const double cw_coef = val_array[T.ElementNo];
    double x[3];
    mfem::Vector transip(x, 3);
    T.Transform(ip, transip);
    const double func_val = (*Function)(transip, param);
    return cw_coef * func_val;
  }

protected:
  double(*Function)(const mfem::Vector&, const Parameters&);
  const Parameters& param;
};



class RhoCoefficient : public mfem::Coefficient
{
public:
  RhoCoefficient(const Parameters& p)
    : param(p)
  { }
  virtual ~RhoCoefficient() { }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip);

protected:
  const Parameters& param;
};

class RhoFuncCoefficient : public mfem::Coefficient
{
public:
  RhoFuncCoefficient(double(*F)(const mfem::Vector&, const Parameters&),
                     const Parameters& p)
    : Function(F), param(p)
  { }
  virtual ~RhoFuncCoefficient() { }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip);

protected:
  double(*Function)(const mfem::Vector&, const Parameters&);
  const Parameters& param;
};

class LambdaFuncCoefficient : public mfem::Coefficient
{
public:
  LambdaFuncCoefficient(double(*F)(const mfem::Vector&, const Parameters&),
                        const Parameters& p)
    : Function(F), param(p)
  { }
  virtual ~LambdaFuncCoefficient() { }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip);

protected:
  double(*Function)(const mfem::Vector&, const Parameters&);
  const Parameters& param;
};

class MuFuncCoefficient : public mfem::Coefficient
{
public:
  MuFuncCoefficient(double(*F)(const mfem::Vector&, const Parameters&),
                    const Parameters& p)
    : Function(F), param(p)
  { }
  virtual ~MuFuncCoefficient() { }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip);

protected:
  double(*Function)(const mfem::Vector&, const Parameters&);
  const Parameters& param;
};



void show_SRM_damp_weights(const Parameters& param);

mfem::Vector compute_function_at_point(double sx, double sy, int nx, int ny,
                                       const mfem::Mesh& mesh,
                                       const mfem::Vertex& point, int cell,
                                       const mfem::GridFunction& U);

mfem::Vector compute_function_at_points(double sx, double sy, int nx, int ny,
                                        const mfem::Mesh& mesh,
                                        const std::vector<mfem::Vertex>& points,
                                        const std::vector<int>& cells_containing_points,
                                        const mfem::GridFunction& U);

void output_snapshots(int time_step, const std::string& snapshot_filebase,
                      const Parameters& param, const mfem::GridFunction& U,
                      const mfem::GridFunction& V);

void output_seismograms(const Parameters& param, const mfem::Mesh& mesh,
                        const mfem::GridFunction &U, const mfem::GridFunction &V,
                        std::ofstream *seisU, std::ofstream *seisV);

#endif // ELWAVE2D_HPP
