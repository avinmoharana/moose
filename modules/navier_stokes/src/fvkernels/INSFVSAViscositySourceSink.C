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
  params.addRequiredParam<MooseFunctorName>(NS::density, "density");
  params.addRequiredParam<MooseFunctorName>("Cv1", "Coef.");
  params.addRequiredParam<MooseFunctorName>("kappa", "Coef.");
  params.addRequiredParam<MooseFunctorName>("Cb1", "Coef.");
  params.addRequiredParam<MooseFunctorName>("Cb2", "Coef.");
  params.addRequiredParam<MooseFunctorName>("sigma_nu", "Coef.");
  params.addRequiredParam<MooseFunctorName>("Cw3", "Coef.");
  params.addRequiredParam<MooseFunctorName>("Cw2", "Coef.");
  params.set<unsigned short>("ghost_layers") = 3;
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
    _rho(getFunctor<ADReal>(NS::density)),
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

ADReal INSFVSAViscositySourceSink::getSAStrainTensorNormDeformation()
{
  constexpr Real offset = 1e-15; // prevents explosion of sqrt(x) derivative to infinity
  const auto xi = _var(makeElemArg(_current_elem),determineState())
	  / (_mu(makeElemArg(_current_elem),determineState()) / _rho(makeElemArg(_current_elem),determineState()));
  const auto xi_cube = std::pow(xi,3);
  const auto cv1_cube = std::pow(_C_v1(makeElemArg(_current_elem),determineState()),3);
  const auto fv1 = xi_cube/(xi_cube+cv1_cube);
  const auto fv2 = 1. - xi/(1. + xi * fv1);

  const auto & grad_u = _u_var->adGradSln(_current_elem, determineState(), true);
  ADReal symmetric_strain_tensor_norm = 0.;
  ADReal vorticity = 0.;
  ADReal meanStrainRate = Utility::pow<2>(grad_u(0));
  if (_dim >= 2)
  {
    const auto & grad_v = _v_var->adGradSln(_current_elem, determineState(), true);
  
    auto S12 = 0.5*(grad_u(1) + grad_v(0));
    auto S21 = S12;
    auto S22 = grad_v(1);
    meanStrainRate +=  (Utility::pow<2>(S12)+ Utility::pow<2>(S21) + Utility::pow<2>(S22));
    auto omega12 = 0.5*(grad_u(1) - grad_v(0));
    auto omega21 = -omega12;
    vorticity +=  (Utility::pow<2>(omega12)+ Utility::pow<2>(omega21) );

    if (_dim >= 3)
    {
      const auto & grad_w = _w_var->adGradSln(_current_elem, determineState(), true);
      auto S13 = 0.5*(grad_u(2) + grad_w(0));
      auto S31 = S13;
      auto S23 = 0.5*(grad_u(2) + grad_w(1));
      auto S32 = S23;
      auto S33 = grad_w(2);

      meanStrainRate +=  (Utility::pow<2>(S13)+ Utility::pow<2>(S31) + Utility::pow<2>(S33));
      meanStrainRate +=  (Utility::pow<2>(S23)+ Utility::pow<2>(S32));
      auto omega13 = 0.5*(grad_u(2) - grad_w(0));
      auto omega23 = 0.5*(grad_u(2) - grad_w(1));
      auto omega31 = -omega13;
      auto omega32 = -omega23;
      
      vorticity +=  (Utility::pow<2>(omega13)+ Utility::pow<2>(omega31) );
      vorticity +=  (Utility::pow<2>(omega23)+ Utility::pow<2>(omega32) );
    }
  }

  ADReal sqrtVorticity = std::sqrt(2.0*vorticity+offset);
  ADReal sqrtStrainRate = std::sqrt(2.0*meanStrainRate+offset);

  const auto kd_sq =std::pow( _kappa(makeElemArg(_current_elem),determineState())*_d(makeElemArg(_current_elem),determineState()), 2);
  const auto S_add = _var(makeElemArg(_current_elem),determineState()) *fv2 / kd_sq;

  // From ANSYS Fluent Theory Manual
  ADReal S_tilda = sqrtVorticity + S_add;

  // Alternate formulation ANSYS Fluent Theory Manual
  //ADReal S_tilda = sqrtVorticity + S_add + 2.0*std::min(0.,sqrtStrainRate-sqrtVorticity);
  //ADReal S_tilda = sqrtVorticity + S_add + 1.0*std::min(0.,sqrtStrainRate-sqrtVorticity);

  // Limiting S_tilda
  //S_tilda = (S_tilda < 0.3*sqrtStrainRate) ? 0.3*sqrtStrainRate : S_tilda; // Min Shat is 0.3 times vorticity
  S_tilda = (S_tilda < 0.3*sqrtVorticity) ? 0.3*sqrtVorticity : S_tilda; // Min Shat is 0.3 times vorticity

  //return(S_tilda+1e-6);
  return(S_tilda);
}

ADReal INSFVSAViscositySourceSink::production()
{
  constexpr Real protection_nu_bar = 0.0000;
  const auto S_tilda = getSAStrainTensorNormDeformation();
  ADReal prod = _C_b1(makeElemArg(_current_elem),determineState())
	  * S_tilda 
	  * _var(makeElemArg(_current_elem),determineState()) + protection_nu_bar;
  return (prod);
}

ADReal INSFVSAViscositySourceSink::destruction()
{
  constexpr Real protection_nu_bar = 0.;

  auto S_tilda = getSAStrainTensorNormDeformation();
  S_tilda = (S_tilda < 1e-6) ? 1e-6 : S_tilda;
  const auto kd_sq =std::pow( _kappa(makeElemArg(_current_elem),determineState())*_d(makeElemArg(_current_elem),determineState()), 2);
  
  ADReal r =  _var(makeElemArg(_current_elem),determineState()) / S_tilda / kd_sq;
  //r = (r > 0.0) ? r : 1e-3;
  r = (r < 10.0) ? r : 10.0;  // limiting r to 10, literature reports r < 1
  const auto g = r + _C_w2(makeElemArg(_current_elem),determineState()) * (std::pow(r,6) - r);
  const auto cw3_6 = std::pow(_C_w3(makeElemArg(_current_elem),determineState()),6);
  const auto fw = g * std::pow( ( 1.0 + cw3_6) / (std::pow(g,6)+cw3_6),1.0/6.0);
  const auto Cw1 = _C_b1(makeElemArg(_current_elem),determineState()) / std::pow(_kappa(makeElemArg(_current_elem),determineState()),2)
             + ( 1.0 + _C_b2(makeElemArg(_current_elem),determineState()) ) / _sigma_nu(makeElemArg(_current_elem),determineState());
  const auto nubar_by_d_sq =std::pow( _var(makeElemArg(_current_elem),determineState())/_d(makeElemArg(_current_elem),determineState()),2);

  ADReal destruct = Cw1 * fw *nubar_by_d_sq + protection_nu_bar;

  return (destruct);
}

ADReal INSFVSAViscositySourceSink::gradSquareTerm()
{
  auto dVar = _var.gradient(makeElemArg(_current_elem),determineState());
  ADReal gst =  _C_b2(makeElemArg(_current_elem),determineState())
	 * (  std::pow(dVar(0),2) + std::pow(dVar(1),2) )
  	 / _sigma_nu(makeElemArg(_current_elem),determineState());
  
  return gst;
}

ADReal
INSFVSAViscositySourceSink::computeQpResidual()
{
#ifdef MOOSE_GLOBAL_AD_INDEXING

  ADReal residual = 0.0;

  ADReal prod = production();
  ADReal destruct = destruction();
  ADReal gst = gradSquareTerm();

  ADReal pt = prod + gst;
  //positiveTerms = positiveTerms > 0 ? positiveTerms : 0.0;
  //destruct = destruct > 0 ? destruct : 0.0;

  //ADReal n = 5.0;
  //if (std::abs(pt)/std::abs(destruct) > n)
  //  pt = std::pow(n*std::abs(destruct.value())+1e-10,2./3.);
  
  //residual += pt - destruct;
  residual += destruct - pt;

  ADReal one = 1.0;

  //return 0;
  
  return residual;

  #else
    return 0;

  #endif
}
