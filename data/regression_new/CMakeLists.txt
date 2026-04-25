# Regression tests for testing analytics library
# Layout: analytics/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/analytics")

# Regression tests for testing boundary_conditions
# Layout: boundary_conditions/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions")

# Regression tests for testing constraints
# Layout: boundary_conditions/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/constraints/generic_vec3")

# Regression tests for testing constraints
# Layout: boundary_conditions/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/constraints/walls")

# Regression tests for testing domain-related operators
# Layout: boundary_conditions/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/domain")

# Regression tests for testing local average
# Layout: compute_local_average/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/compute_local_average")

# Regression tests for testing compute local strain
# Layout: compute_local_strain/<variant>.msp
# Comment for now : this depends on onika bugfix-xs_data_move branch
#AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/compute_local_strain")

# Regression tests for numerical schemes (nve, nvt, npt)
# Layout: numerical_schemes/<testcase>/<variant>.msp
AddRegressionTestSimple("${CMAKE_CURRENT_SOURCE_DIR}/numerical_schemes")

# Regression tests for flexible molecules
# Layout: flexible_molecule/<testcase>/<variant>.msp
AddRegressionTestSimple("${CMAKE_CURRENT_SOURCE_DIR}/flexible_molecule")

# Regression tests for pair potentials
# Layout: potentials/pair/<testcase>/<variant>.msp
AddRegressionTestSimple("${CMAKE_CURRENT_SOURCE_DIR}/potentials/pair")

# Regression tests for eam potentials
# Layout: potentials/eam/<testcase>/<variant>.msp
AddRegressionTestSimple("${CMAKE_CURRENT_SOURCE_DIR}/potentials/eam")

# AddRegressionTestDir("${CMAKE_SOURCE_DIR}/data/regression_new/potentials/pair")
# AddRegressionTestDir("${CMAKE_SOURCE_DIR}/data/regression_new/potentials/eam")

# # Regression test directory for numerical schemes (nve, nvt, npt)
# AddRegressionTestSimple("${CMAKE_CURRENT_SOURCE_DIR}/numerical_schemes")
# # Regression test directory for flexible molecules (nve, nvt, npt)
# AddRegressionTestSimple("${CMAKE_CURRENT_SOURCE_DIR}/flexible_molecule")

# #  AddRegressionTestDir("${CMAKE_SOURCE_DIR}/data/regression_new/potentials/pair")
# #  AddRegressionTestDir("${CMAKE_SOURCE_DIR}/data/regression_new/potentials/eam")
