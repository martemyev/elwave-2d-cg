#ifndef ELWAVE2D_HPP
#define ELWAVE2D_HPP

#include "config.hpp"
#include "mfem.hpp"

#include <vector>

class Parameters;
class Source;



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



/**
 * Implementation of a vector point force as a component of a source term.
 */
class VectorPointForce: public mfem::VectorCoefficient
{
public:
  VectorPointForce(int dim, const Source& s);
  ~VectorPointForce() { }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip);

private:
  const Source& source;
};



/**
 * Implementation of a moment tensor density as a component of a source term.
 */
class MomentTensorSource: public mfem::VectorCoefficient
{
public:
  MomentTensorSource(int dim, const Source& s);
  ~MomentTensorSource() { }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip);

private:
  const Source& source;
};




class ElasticWave2D
{
public:
  ElasticWave2D(const Parameters& _param);
  ~ElasticWave2D();

  void run();

private:
  const Parameters& param;

  /*mfem::Mesh *mesh;
  mfem::FiniteElementCollection *fec;
  mfem::FiniteElementSpace *fespace;
  mfem::BilinearForm *stif;
  mfem::BilinearForm *mass;
  mfem::BilinearForm *dampM, *dampS;

  CWConstCoefficient *rho_coef;
  CWConstCoefficient *lambda_coef;
  CWConstCoefficient *mu_coef;
  CWConstCoefficient *rho_damp_coef;
  CWConstCoefficient *lambda_damp_coef;
  CWConstCoefficient *mu_damp_coef;

  mfem::ElasticityIntegrator *elast_int;
  mfem::VectorMassIntegrator *mass_int;
  mfem::VectorMassIntegrator *damp_int;
  VectorPointForce *vector_point_force;
  MomentTensorSource *momemt_tensor_source;
  mfem::VectorDomainLFIntegrator *point_force_int;
  mfem::VectorDomainLFIntegrator *moment_tensor_int;

  mfem::LinearForm *b;*/

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



mfem::Vector compute_solution_at_points(const std::vector<mfem::Vertex>& points,
                                        const std::vector<int>& cells_containing_points,
                                        const mfem::GridFunction& U);

void show_SRM_damp_weights(const Parameters& param);

#endif // ELWAVE2D_HPP
