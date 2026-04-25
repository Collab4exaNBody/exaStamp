# Regression tests for testing analytics library
# Layout: analytics/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/analytics")

# Regression tests for testing boundary_conditions
# Layout: boundary_conditions/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions")

# Regression tests for testing generic_vec3 operators
# Layout: boundary_conditions/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/constraints/generic_vec3")

# Regression tests for testing wall-related operators
# Layout: boundary_conditions/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/constraints/walls")

# Regression tests for testing deterministic noise related operators
# Layout: boundary_conditions/<variant>.msp
AddRegressionTestSimple("${CMAKE_CURRENT_SOURCE_DIR}/deterministic")

# Regression tests for testing domain-related operators
# Layout: boundary_conditions/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/domain")

# Regression tests for testing local average
# Layout: compute_local_average/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/compute_local_average")

# Regression tests for testing compute local strain
# Layout: compute_local_strain/<variant>.msp
# WARNING : this depends on onika bugfix-xs_data_move branch
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/compute_local_strain")



# Regression tests for flexible molecules
# Layout: flexible_molecule/<testcase>/<variant>.msp
#AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/flexible_molecule/molbonds")
#AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/flexible_molecule/molperf")
#AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/flexible_molecule/molstampv4")
#AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/flexible_molecule/molxyz")
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/flexible_molecule/polybutadiene")
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/flexible_molecule/tetrahydrofurane")
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/flexible_molecule/triaminotrinitrobenzene")

# Regression tests for testing lattice operators
# Layout: compute_local_strain/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/lattice")


# Regression tests for numerical schemes (nve, nvt, npt)
# Layout: numerical_schemes/<testcase>/<variant>.msp
AddRegressionTestSimple("${CMAKE_CURRENT_SOURCE_DIR}/numerical_schemes")

# Regression tests for testing particle regions
# Layout: compute_local_strain/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/particle_regions")

# Regression tests for pair potentials
# Layout: potentials/pair/<testcase>/<variant>.msp
AddRegressionTestSimple("${CMAKE_CURRENT_SOURCE_DIR}/potentials/pair")

# Regression tests for eam potentials
# Layout: potentials/eam/<testcase>/<variant>.msp
AddRegressionTestSimple("${CMAKE_CURRENT_SOURCE_DIR}/potentials/eam")

# Regression tests for eam potentials
# Layout: potentials/eam/<testcase>/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/random")

# Regression tests for applied deformation using xform_function
# Layout: potentials/eam/<testcase>/<variant>.msp
AddRegressionTestFlat("${CMAKE_CURRENT_SOURCE_DIR}/xform_function")