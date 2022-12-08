//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVSAViscositySourceSink.h"
#include "NS.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", INSFVSAViscositySourceSink);

InputParameters
INSFVSAViscositySourceSink::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription(
      "Elemental kernel to compute the production and destruction "
      " terms of Spalart-Allmaras turbulence model.");
  params.addRequiredCoupledVar("u", "The velocity in the x direction.");
  params.addCoupledVar("v", "The velocity in the y direction.");
  params.addCoupledVar("w", "The velocity in the z direction.");
  params.addParam<std::vector<BoundaryName>>("walls", "Boundaries corresponding to solid walls");
  params.addRequiredParam<MooseFunctorName>("distance", "distance");
  params.addRequiredParam<MooseFunctorName>("mu", "viscosity");
  params.addRequiredParam<MooseFunctorName>("rho", "density");
  params.addRequiredParam<MooseFunctorName>("Cv1", "Coef.");
  params.addRequiredParam<MooseFunctorName>("kappa", "Coef.");
  params.addRequiredParam<MooseFunctorName>("Cb1", "Coef.");
  params.addRequiredParam<MooseFunctorName>("Cb2", "Coef.");
  params.addRequiredParam<MooseFunctorName>("sigma_nu", "Coef.");
  params.addRequiredParam<MooseFunctorName>("Cw3", "Coef.");
  params.addRequiredParam<MooseFunctorName>("Cw2", "Coef.");
  params.set<unsigned short>("ghost_layers") = 2;
  return params;
}

INSFVSAViscositySourceSink::INSFVSAViscositySourceSink(const InputParameters & params)
  : FVElementalKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _u_var(dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("u", 0))),
    _v_var(params.isParamValid("v")
               ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("v", 0))
               : nullptr),
    _w_var(params.isParamValid("w")
               ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("w", 0))
               : nullptr),
    _d(getFunctor<ADReal>("distance")),	       
    _mu(getFunctor<ADReal>("mu")),
    _rho(getFunctor<ADReal>("rho")),
    _C_v1(getFunctor<ADReal>("Cv1")),
    _kappa(getFunctor<ADReal>("kappa")),
    _C_b1(getFunctor<ADReal>("Cb1")),
    _C_b2(getFunctor<ADReal>("Cb2")),
    _sigma_nu(getFunctor<ADReal>("sigma_nu")),
    _C_w3(getFunctor<ADReal>("Cw3")),
    _C_w2(getFunctor<ADReal>("Cw2"))
{
#ifndef MOOSE_GLOBAL_AD_INDEXING
  mooseError("INSFV is not supported by local AD indexing. In order to use INSFV, please run the "
             "configure script in the root MOOSE directory with the configure option "
             "'--with-ad-indexing-type=global'");
#endif

  if (!_u_var)
    paramError("u", "the u velocity must be an INSFVVelocityVariable.");

  if (_dim >= 2 && !_v_var)
    paramError("v",
               "In two or more dimensions, the v velocity must be supplied and it must be an "
               "INSFVVelocityVariable.");

  if (_dim >= 3 && !_w_var)
    paramError("w",
               "In three-dimensions, the w velocity must be supplied and it must be an "
               "INSFVVelocityVariable.");
}
ADReal INSFVSAViscositySourceSink::getSymmetricStrainTensorNorm()
{
  constexpr Real offset = 1e-15; // prevents explosion of sqrt(x) derivative to infinity

  const auto & grad_u = _u_var->adGradSln(_current_elem);
  ADReal symmetric_strain_tensor_norm = 0.;
  ADReal deformation = 2.0 * Utility::pow<2>(grad_u(0));
  ADReal vorticity = 0.;
  if (_dim >= 2)
  {
    const auto & grad_v = _v_var->adGradSln(_current_elem);
    deformation +=
        2.0*Utility::pow<2>(grad_v(1)) + Utility::pow<2>(grad_v(0) + grad_u(1));
    vorticity += 2.0 * Utility::pow<2>(grad_v(0) - grad_u(1));

    if (_dim >= 3)
    {
      const auto & grad_w = _w_var->adGradSln(_current_elem);
      deformation += 2.0 * Utility::pow<2>(grad_w(2)) +
                     Utility::pow<2>(grad_u(2) + grad_w(0)) +
		     Utility::pow<2>(grad_v(2) + grad_w(1));
      // TODO vorticity 3D
    }
  }

  return(std::sqrt(vorticity+offset) 
         + 2.0 * std::min(0.,std::sqrt(deformation+offset)-std::sqrt(vorticity+offset)));

}

ADReal INSFVSAViscositySourceSink::production()
{
  constexpr Real protection_nu_bar = 0.001;
  
  auto omega = getSymmetricStrainTensorNorm();
  auto xi = _var(makeElemArg(_current_elem))
	  / (_mu(makeElemArg(_current_elem)) / _rho(makeElemArg(_current_elem)));
  auto fv1 = Utility::pow<3>(xi) 
	  / ( Utility::pow<3>(xi) + Utility::pow<3>(_C_v1(makeElemArg(_current_elem))) );
  auto fv2 = 1. - xi/(1. + xi * fv1);
  auto kd_sq = Utility::pow<2>(_kappa(makeElemArg(_current_elem))
              *_d(makeElemArg(_current_elem)));
  auto S_tilda = omega + _var(makeElemArg(_current_elem)) *fv2 / kd_sq;
  S_tilda = std::max(S_tilda, 0.3*omega);
  auto prod = (_C_b1(makeElemArg(_current_elem)) 
	  * S_tilda 
	  * _var(makeElemArg(_current_elem))) + protection_nu_bar;
  
  prod = (prod > 0.) ? prod : 0.0;
  return (prod);
}

ADReal INSFVSAViscositySourceSink::destruction()
{
  constexpr Real protection_nu_bar = 0.001;

  auto omega = getSymmetricStrainTensorNorm();

  auto kd_sq = Utility::pow<2>(_kappa(makeElemArg(_current_elem))
                               *_d(makeElemArg(_current_elem)));
  auto xi = _var(makeElemArg(_current_elem))
	  / ( _mu(makeElemArg(_current_elem)) / _rho(makeElemArg(_current_elem)) );
  auto fv1 = Utility::pow<3>(xi) 
	  / ( Utility::pow<3>(xi) + Utility::pow<3>(_C_v1(makeElemArg(_current_elem))) );
  auto fv2 = 1. - xi/(1. + xi * fv1);
  auto S_tilda = omega + _var(makeElemArg(_current_elem)) *fv2 / kd_sq;
  S_tilda = std::max(S_tilda, 0.3*omega);   // OpenFOAM implementation, ensures minimum positive value
  auto r =  _var(makeElemArg(_current_elem)) / S_tilda / kd_sq;
  r = (r < 10.0) ? r : 10.0;  // OpenFOAM implementation, literature says r < 1
  auto g = r + _C_w2(makeElemArg(_current_elem)) * (Utility::pow<6>(r) - r);
  auto cw3_6 = Utility::pow<6>(_C_w3(makeElemArg(_current_elem)));
  auto fw = g * std::pow( ( 1.0 + cw3_6) / (Utility::pow<6>(g)+cw3_6),1.0/6.0);
  auto Cw1 = _C_b1(makeElemArg(_current_elem)) / Utility::pow<2>(_kappa(makeElemArg(_current_elem)))
             + ( 1.0 + _C_b2(makeElemArg(_current_elem)) ) / _sigma_nu(makeElemArg(_current_elem));
  auto nubar_by_d_sq = Utility::pow<2>( _var(makeElemArg(_current_elem))/_d(makeElemArg(_current_elem)) );	     
  ADReal destruct = Cw1 * fw * nubar_by_d_sq + protection_nu_bar;
  
  destruct = (destruct > 0.) ? destruct : 0.0;
  return (destruct);
}

ADReal INSFVSAViscositySourceSink::gradSquareTerm()
{
  auto dVar = _var.gradient(makeElemArg(_current_elem));
  return ( _C_b2(makeElemArg(_current_elem))
	 / _sigma_nu(makeElemArg(_current_elem))
	 * (  std::pow(dVar(0),2) + std::pow(dVar(1),2) + std::pow(dVar(2),2)));
}

ADReal
INSFVSAViscositySourceSink::computeQpResidual()
{
#ifdef MOOSE_GLOBAL_AD_INDEXING

  ADReal residual = 0.0;

  ADReal prod = production();
  ADReal destruct = destruction();
  ADReal gst = gradSquareTerm(); 
  
  residual += -prod + destruct - gst;

  return residual;

  #else
    return 0;

  #endif
}
