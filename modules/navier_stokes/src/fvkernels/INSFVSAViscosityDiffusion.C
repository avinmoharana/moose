//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVSAViscosityDiffusion.h"
#include "INSFVEnergyVariable.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", INSFVSAViscosityDiffusion);

InputParameters
INSFVSAViscosityDiffusion::validParams()
{
  auto params = FVFluxKernel::validParams();
  params.addClassDescription("Diffusion term in the SA turbulence model"
                             " $-1/sigma_nu * div(mu_t/rho * grad(nu))$");

  params.addRequiredParam<MooseFunctorName>("mu_t", "turbulent viscosity");
  params.addRequiredParam<MooseFunctorName>("rho", "density");
  params.addRequiredParam<MooseFunctorName>("sigma_nu", "coefficient in diffusion term");
  
  params.set<unsigned short>("ghost_layers") = 2;
  return params;
}

INSFVSAViscosityDiffusion::INSFVSAViscosityDiffusion(const InputParameters & params)
  : FVFluxKernel(params),
    _mu_t(getFunctor<ADReal>("mu_t")),
    _rho(getFunctor<ADReal>("rho")),
    _sigma_nu(getFunctor<ADReal>("sigma_nu"))
{
#ifndef MOOSE_GLOBAL_AD_INDEXING
  mooseError("INSFV is not supported by local AD indexing. In order to use INSFV, please run "
             "the configure script in the root MOOSE directory with the configure option "
             "'--with-ad-indexing-type=global'");
#endif
  if (!dynamic_cast<INSFVEnergyVariable *>(&_var))
    mooseError("INSFVSAViscosityDiffusion may only be used with a fluid temperature variable, "
               "of variable type INSFVEnergyVariable.");
}

ADReal
INSFVSAViscosityDiffusion::computeQpResidual()
{
  // Interpolate effective diffusity on the face
  ADReal diffusion_face;
  const auto face_elem = elemArg();
  const auto face_neighbor = neighborArg();

  const auto mut_fe = _mu_t(face_elem,determineState());
  const auto rho_fe = _rho(face_elem,determineState());
  const auto sigma_fe = _sigma_nu(face_elem,determineState());
  const auto mut_fn = _mu_t(face_neighbor,determineState());
  const auto rho_fn = _rho(face_neighbor,determineState());
  const auto sigma_fn = _sigma_nu(face_neighbor,determineState());
 
 Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
		  	diffusion_face, mut_fe/rho_fe*sigma_fe, mut_fn/rho_fn*sigma_fn, *_face_info, true);
  
 //std::cout<<"cp2"<<std::endl;
  // Compute nu gradient dotted with the surface normal
  auto dVardn = gradUDotNormal(determineState());

 //std::cout<<"cp3"<<std::endl;
  return -diffusion_face * dVardn;
}
