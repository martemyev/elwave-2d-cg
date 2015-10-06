#ifndef ELWAVE2D_HPP
#define ELWAVE2D_HPP

#include "mfem.hpp"

class Parameters;
class Source;



/**
 * Cell-wise constant coefficient
 */
class CWConstCoefficient : public mfem::Coefficient
{
public:
  CWConstCoefficient(double *array, bool own = 1)
    : val_array(array), own_array(own)
  { }

  ~CWConstCoefficient() { if (own_array) delete[] val_array; }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
  {
    return val_array[T.ElementNo];
  }

private:
  double *val_array;
  bool own_array;
};



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

  enum Method { FEM, SEM };

  void run_FEM();
  void run_SEM();


private:
  const Parameters& param;

  mfem::Mesh *mesh;
  mfem::FiniteElementCollection *fec;
  mfem::FiniteElementSpace *fespace;
  mfem::BilinearForm *stif;
  mfem::BilinearForm *mass;

  mfem::ElasticityIntegrator *elast_int;
  mfem::VectorMassIntegrator *mass_int;
  VectorPointForce *vector_point_force;
  MomentTensorSource *momemt_tensor_source;
  mfem::VectorDomainLFIntegrator *point_force_int;
  mfem::VectorDomainLFIntegrator *moment_tensor_int;

  mfem::LinearForm *b;

  void offline_stage();
  void online_stage(Method method);
};

#endif // ELWAVE2D_HPP
