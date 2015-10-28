#ifndef SOURCE_HPP
#define SOURCE_HPP

#include "config.hpp"
#include "mfem.hpp"



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
  double GaussFirstDerivative(double t) const;
  void PointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void MomentTensorSource(const mfem::Vector& x, mfem::Vector& f) const;

private:
  void DeltaPointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void GaussPointForce(const mfem::Vector& x, mfem::Vector& f) const;
  void DivDeltaMomentTensor(const mfem::Vector& x, mfem::Vector& f) const;
  void DivGaussMomentTensor(const mfem::Vector& x, mfem::Vector& f) const;
};

#endif // SOURCE_HPP
