//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

class INSFVVelocityVariable;

/*
 * Computes the value of the eddy viscosity for the mixing length model.
 */
class spalartAllmarasViscosity : public AuxKernel
{
public:
  static InputParameters validParams();

  spalartAllmarasViscosity(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Spalart-Allmaras equivalent kinematic viscosity
  const Moose::Functor<ADReal> & _nu;

  /// Fluid Density
  const Moose::Functor<ADReal> & _rho;

  /// Dynamic Viscosity
  const Moose::Functor<ADReal> & _mu;

  /// C-v1 closure constant
  const Moose::Functor<ADReal> & _C_v1;

};
