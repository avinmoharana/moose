mu = 0.0000181
rho = 1.225

bulk_u = 20.0

# S-A standard coefficients
C_b1 = 0.1355
C_b2 = 0.622
sigma_nu = 0.6667
C_v1 = 7.1
C_w2 = 0.3
C_w3 = 2.0
kappa = 0.4187

advected_interp_method = 'upwind'
velocity_interp_method = 'rc'

[Mesh]
  file = 2dTankMesh.e
[]

[Outputs]
  exodus = true
  [debug]
    type = Exodus
    execute_on = nonlinear
  []
[]

[Problem]
  fv_bcs_integrity_check = true
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  # The upwind and Rhie-Chow interpolation schemes are used here.
  # advected_interp_method='upwind'
  # velocity_interp_method='rc'
  two_term_boundary_expansion = false
  kernel_coverage_check = false
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = u
    v = v
    pressure = pressure
  []
[]

[Variables]
  [u]
    type = INSFVVelocityVariable
    initial_condition = 0.00001
  []
  [v]
    type = INSFVVelocityVariable
    initial_condition = 1e-6
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [nu_bar]
    type = INSFVEnergyVariable
    initial_condition = 1
  []
[]


[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = u
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = u
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_viscosity_rans]
    type = INSFVMomentumDiffusion
    variable = u
    mu = mu_t
    momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = u
    momentum_component = 'x'
    pressure = pressure
  []

  [v_advection]
    type = INSFVMomentumAdvection
    variable = v
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = v
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_viscosity_rans]
    type = INSFVMomentumDiffusion
    variable = v
    mu = mu_t
    momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = v
    momentum_component = 'y'
    pressure = pressure
  []
  [SA_source_sink]
    type = INSFVSAViscositySourceSink
    variable = nu_bar
    u = u
    v = v
    distance = d
    mu = ${mu}
    rho = ${rho}
    Cv1 = ${C_v1}
    kappa = ${kappa}
    Cb1 = ${C_b1}
    Cb2 = ${C_b2}
    sigma_nu = ${sigma_nu}
    Cw3 = ${C_w3}
    Cw2 = ${C_w2}
  []
  [SA_diffusion]
    type = INSFVSAViscosityDiffusion
    variable = nu_bar
    mu_t = mu_t
    rho = ${rho}
    sigma_nu = ${sigma_nu}
  []
  [SA_diffusion_laminar]
    type = INSFVSAViscosityDiffusion
    variable = nu_bar
    mu_t = ${mu}
    rho = ${rho}
    sigma_nu = ${sigma_nu}
  []
  [SA_advection]
    type = INSFVSAViscosityAdvection
    variable = nu_bar
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
  []

[]

[AuxVariables]
  [mu_t]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    initial_condition = 1
  []
  [d]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [compute_mu_t]
    type = spalartAllmarasViscosity
    variable = mu_t
    nu = nu_bar
    rho = ${rho}
    mu = ${mu}
    C_v1 = ${C_v1}
  []
  [compute_distance]
    type = WallDistanceAux
    walls = 'WallSurface'
    variable = d
    execute_on = 'initial'
  []
[]

[FVBCs]
  [inlet_u]
    type = INSFVInletVelocityBC
    boundary = 'Orifice'
    variable = u
    function = ${bulk_u}
    #function = '20*t'
  []
  [inlet_v]
    type = INSFVInletVelocityBC
    boundary = 'Orifice'
    variable = v
    function = '0'
  []
  [inlet_nu]
    type = FVDirichletBC
    boundary = 'Orifice'
    variable = nu_bar
    value = 0.001
  []
  [walls_u]
    type = INSFVNoSlipWallBC
    boundary = 'WallSurface'
    variable = u
    function = 0
  []
  [walls_v]
    type = INSFVNoSlipWallBC
    boundary = 'WallSurface'
    variable = v
    function = 0
  []
  [walls_nu]
    type = INSFVNoSlipWallBC
    boundary = 'WallSurface'
    variable = nu_bar
    function = 0
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'Outlet'
    variable = pressure
    function = 0
  []
[]
[Preconditioning]
  [SMP]
    type = SMP
    full = false
  []
[]

[Executioner]

 type = Steady
 solve_type = 'PJFNK'
 petsc_options_iname = '-pc_type'
 petsc_options_value = 'lu'
 line_search = 'none'
 #automatic_scaling = true
 nl_rel_tol = 1e-3
 nl_abs_tol = 1e-5

 #type = Steady
 #solve_type = 'NEWTON'
 #petsc_options_iname = '-pc_type'
 #petsc_options_value = 'lu '
 #line_search = 'none'
 #nl_rel_tol = 1e-6

 #type = Transient
 #solve_type = 'PJFNK'
 #petsc_options_iname = '-pc_type'
 #petsc_options_value = 'lu'
 #automatic_scaling = true
 #line_search = 'none'
 #[TimeStepper]
 #  type = IterationAdaptiveDT
 #  optimal_iterations = 10
 #  dt = 1e-2
 #[]
 #nl_abs_tol = 1e-8
 #end_time = 1
[]

