//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
 * SA source term implementation follows this reference
 * https://turbmodels.larc.nasa.gov/spalart.html
*/

#pragma once

#include "FVElementalKernel.h"

/**
 * Simple class to demonstrate off diagonal Jacobian contributions.
 */
class INSFVSAViscositySourceSink : public FVElementalKernel
{
public:
  static InputParameters validParams();

  INSFVSAViscositySourceSink(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

protected:
  virtual ADReal getSAStrainTensorNormDeformation();
  virtual ADReal getSAStrainTensorNormVorticity();
  virtual ADReal production();
  virtual ADReal destruction();
  virtual ADReal gradSquareTerm();

  /// The dimension of the simulation
  const unsigned int _dim;

  /// x-velocity
  const INSFVVelocityVariable * const _u_var;
  /// y-velocity
  const INSFVVelocityVariable * const _v_var;
  /// z-velocity
  const INSFVVelocityVariable * const _w_var;
  /// distance
  const Moose::Functor<ADReal> &_d;
  /// visosity
  const Moose::Functor<ADReal> &_mu;
  /// Density
  const Moose::Functor<ADReal> &_rho;

  /// Cv1
  const Moose::Functor<ADReal> & _C_v1;
  /// kappa
  const Moose::Functor<ADReal> & _kappa;
  /// Cb1
  const Moose::Functor<ADReal> & _C_b1;
  /// Cb2
  const Moose::Functor<ADReal> & _C_b2;
  /// sigma_nu
  const Moose::Functor<ADReal> & _sigma_nu;
  /// Cw3
  const Moose::Functor<ADReal> & _C_w3;
  /// Cw2
  const Moose::Functor<ADReal> & _C_w2;

};
