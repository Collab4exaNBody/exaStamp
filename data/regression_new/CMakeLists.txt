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
