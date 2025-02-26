//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GapFluxModelRadiation.h"
#include "libmesh/utility.h"

registerMooseObject("HeatConductionApp", GapFluxModelRadiation);

InputParameters
GapFluxModelRadiation::validParams()
{
  InputParameters params = GapFluxModelBase::validParams();
  params.addClassDescription("Gap flux model with a constant conductance");
  params.addCoupledVar("T", "Temperature");
  params.addParam<Real>("stefan_boltzmann", 5.670373e-8, "Stefan-Boltzmann constant");
  params.addRangeCheckedParam<Real>("primary_emissivity",
                                    1,
                                    "primary_emissivity>=0 & primary_emissivity<=1",
                                    "The emissivity of the primary surface");
  params.addRangeCheckedParam<Real>("secondary_emissivity",
                                    1,
                                    "secondary_emissivity>=0 & secondary_emissivity<=1",
                                    "The emissivity of the secondary surface");
  return params;
}

GapFluxModelRadiation::GapFluxModelRadiation(const InputParameters & parameters)
  : GapFluxModelBase(parameters),
    _primary_T(adCoupledNeighborValue("T")),
    _secondary_T(adCoupledValue("T")),
    _stefan_boltzmann(getParam<Real>("stefan_boltzmann"))
{
  const auto emissivity_primary = getParam<Real>("primary_emissivity");
  const auto emissivity_secondary = getParam<Real>("secondary_emissivity");

  _emissivity = emissivity_primary != 0.0 && emissivity_secondary != 0.0
                    ? 1.0 / emissivity_primary + 1.0 / emissivity_secondary - 1
                    : 0.0;
}

ADReal
GapFluxModelRadiation::computeFlux() const
{
  /*
   Gap conductance due to radiation is based on the diffusion approximation:

      qr = sigma*Fe*(Tf^4 - Tc^4) ~ hr(Tf - Tc)
         where sigma is the Stefan-Boltzmann constant, Fe is an emissivity function, Tf and Tc
         are the fuel and clad absolute temperatures, respectively, and hr is the radiant gap
         conductance. Solving for hr,

      hr = sigma*Fe*(Tf^4 - Tc^4) / (Tf - Tc)
         which can be factored to give:

      hr = sigma*Fe*(Tf^2 + Tc^2) * (Tf + Tc)

   Approximating the fuel-clad gap as infinite parallel planes, the emissivity function is given by:

      Fe = 1 / (1/ef + 1/ec - 1)
  */

  if (_emissivity == 0.0)
    return 0.0;

  // We add 'surface_integration_factor' to account for the surface integration of the conductance
  // due to radiation.

  const ADReal temp_func =
      (_primary_T[_qp] * _primary_T[_qp] + _secondary_T[_qp] * _secondary_T[_qp]) *
      (_primary_T[_qp] + _secondary_T[_qp]);

  return (_primary_T[_qp] - _secondary_T[_qp]) * _stefan_boltzmann * temp_func / _emissivity *
         _surface_integration_factor;
}
