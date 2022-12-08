//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVFluxKernel.h"

/**
 * A flux kernel for diffusion term of the SA turbulence model,
 * using functor material properties
 */
class INSFVSAViscosityDiffusion : public FVFluxKernel
{
public:
  static InputParameters validParams();
  INSFVSAViscosityDiffusion(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  /// turbulent viscosity
  const Moose::Functor<ADReal> & _mu_t;
  /// density
  const Moose::Functor<ADReal> & _rho;
  /// turbulent coefficient in diffusion term - divides diffusion
  const Moose::Functor<ADReal> & _sigma_nu;
};
