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
#include "MortarInterface.h"

/**
 * Base class for creating new nodally-based mortar auxiliary kernels
 *
 */
template <typename ComputeValueType>
class MortarNodalAuxKernelTempl : public AuxKernelTempl<ComputeValueType>, protected MortarInterface
{
public:
  static InputParameters validParams();

  MortarNodalAuxKernelTempl(const InputParameters & parameters);

  /**
   * Computes the value and stores it in the solution vector
   */
  void compute() override;

protected:
  void precalculateValue() override final;

  using AuxKernelTempl<ComputeValueType>::isNodal;
  using AuxKernelTempl<ComputeValueType>::paramError;
  using AuxKernelTempl<ComputeValueType>::mooseError;
  using AuxKernelTempl<ComputeValueType>::_subproblem;
  using AuxKernelTempl<ComputeValueType>::_assembly;
  using AuxKernelTempl<ComputeValueType>::_current_node;
  using AuxKernelTempl<ComputeValueType>::_var;
  using AuxKernelTempl<ComputeValueType>::computeValue;

  /// Whether we're computing on the displaced mesh
  const bool _displaced;

  /// The base finite element problem
  FEProblemBase & _fe_problem;

  /// The mortar segment volume
  Real _msm_volume;
};

typedef MortarNodalAuxKernelTempl<Real> MortarNodalAuxKernel;
typedef MortarNodalAuxKernelTempl<RealVectorValue> VectorMortarNodalAuxKernel;
