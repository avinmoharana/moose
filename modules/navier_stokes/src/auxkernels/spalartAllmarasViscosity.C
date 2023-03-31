//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "spalartAllmarasViscosity.h"
#include "MooseMesh.h"
#include "INSFVMethods.h"

registerMooseObject("NavierStokesApp", spalartAllmarasViscosity);

InputParameters
spalartAllmarasViscosity::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addClassDescription("Computes the turbulent viscosity according to Spalart-Allmaras model.");
  params.addRequiredParam<MooseFunctorName>("nu", "Splart-Allmaras viscosity");
  params.addRequiredParam<MooseFunctorName>("rho", "fluid density");
  params.addRequiredParam<MooseFunctorName>("mu", "fluid dynamic viscosity");
  params.addRequiredParam<MooseFunctorName>("C_v1", "Closure constant");

  return params;
}

spalartAllmarasViscosity::spalartAllmarasViscosity(const InputParameters & params)
  : AuxKernel(params),
    _nu(getFunctor<ADReal>("nu")),
    _rho(getFunctor<ADReal>("rho")),
    _mu(getFunctor<ADReal>("mu")),
    _C_v1(getFunctor<ADReal>("C_v1"))
{
}

Real
spalartAllmarasViscosity::computeValue()
{
  Real cv1 = _C_v1(makeElemArg(_current_elem)).value();
  Real xi = _nu(makeElemArg(_current_elem)).value()
             /(_mu(makeElemArg(_current_elem))/_rho(makeElemArg(_current_elem))).value();
  Real fv1 = Utility::pow<3>(xi) /
             (Utility::pow<3>(xi) + Utility::pow<3>(cv1));
  Real mut = 0.;
  //if (_nu(makeElemArg(_current_elem)).value() < 0.) mut = 0.;
  //else mut = ( _rho(makeElemArg(_current_elem)) *_nu(makeElemArg(_current_elem))).value() * fv1 ;
  mut = ( _rho(makeElemArg(_current_elem)) *_nu(makeElemArg(_current_elem))).value() * fv1 ;
  return (mut);
}
