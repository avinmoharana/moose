# CSVDiff

`CSVDiff` tests compare `CSV` output(s) for the test to a reference in the specified
`gold_dir` folder. `CSV` output(s) contains time steps, scalar variables and postprocessors.
`CSV` output(s) may also be generated by VectorPostprocessors.
This tester picks up any difference between the reference file and the simulation
output, at any time step included in the output.
More documentation for the `CSVDiff` utility may be found [here](CSVDiff.md).

## Options

Test configuration options are added to the `tests` file.

- `abs_zero`: Sets an absolute tolerance, defaults to 1e-10. Both absolute and relative tolerances must
  be met for a test to pass.

- `rel_err`: Sets a relative tolerance, defaults to 5.5e-6.

- `override_columns`: A list of variable names to customize the `CSVDiff` tolerances

- `override_rel_err`: A list of customized relative error tolerances.

- `override_abs_zero`: A list of customized absolute zero tolerances.

- `comparison_file`: Use supplied custom comparison config file.


Other test commands & restrictions may be found in the [TestHarness documentation](TestHarness.md).

## Example test configuration in the MOOSE test suite

In this example, three `CSVDiff` tests are created to test a particular type of preconditioner.
The data stored in the `CSV` output, the L2 error between the numerical and an analytical solution,
is being compared before and after the steady-state solve.

!listing test/tests/preconditioners/vcp/tests