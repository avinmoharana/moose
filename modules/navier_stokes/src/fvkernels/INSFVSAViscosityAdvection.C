//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSFVSAViscosityAdvection.h"

registerMooseObject("NavierStokesApp", INSFVSAViscosityAdvection);

InputParameters
INSFVSAViscosityAdvection::validParams()
{
  auto params = INSFVAdvectionKernel::validParams();
  params.addClassDescription("Advects SA viscosity.");
  return params;
}

INSFVSAViscosityAdvection::INSFVSAViscosityAdvection(const InputParameters & params)
  : INSFVAdvectionKernel(params)
{
}

ADReal
INSFVSAViscosityAdvection::computeQpResidual()
{
  ADReal adv_quant_interface;

  const auto elem_face = elemArg();
  const auto neighbor_face = neighborArg();

  // Velocity interpolation
  const auto v = _rc_vel_provider.getVelocity(_velocity_interp_method, *_face_info,determineState(), _tid);

  // Interpolation of advected quantity
  Moose::FV::interpolate(_advected_interp_method,
		         adv_quant_interface,
			 _var(elem_face, determineState()),
			 _var(neighbor_face, determineState()),
			 v,
			 *_face_info,
			 true);
  
  return _normal * v * adv_quant_interface;
}
