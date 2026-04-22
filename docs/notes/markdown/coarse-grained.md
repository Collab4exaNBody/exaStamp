# Coarse-grained models

ExaStamp is able to run simulations for coarse-grained particles using either dissipation Particle Dynamics (DPD) ofr isothermal situations or DPD with energy conservation (DPDE) that preserves the total energy as the name suggests.

## <a name="simulation-mode"></a> 1. Simulation mode 
## # 
The choice of a simulation mode is mandatory in ExaStamp. Each enables a different set of functionalities. The mode is selected by adding the `mode` keyword with the desired value. In the rest of the documentation, the options requiring a specific mode are indicated and signaled with the symbols associated with the compatible modes.

| Keyword: `mode`      | Description                                            | Symbol   |
| -------------------- | ------------------------------------------------------ | :------: |
| `single_spec_atom`   | Used for classical MD                                  | *        |
| `single_spec_meso`   | Used for coarse-grained dynamics (DPD, DPDE)           | &dagger; |
| `single_spec_smooth` | Used for Smoothed Dissipative Particle Dynamics (SDPD) | &Dagger; |


## <a name="scheme"></a> 2. Scheme 

The choice of a numerical scheme is mandatory. It tells ExaStamp what dynamics to run. The modes compatible with each scheme are indicated in the following table and links to a more thorough documentation on the dynamics and their integration scheme is given. The schemes available to perform classical MD simulations are listed in the following table.

| Keyword: `scheme`        | Description                                                                                                                         | Modes                                   | Link             |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------  | :-------------------------------------: | :--------------: | 
| `dpd_splitting`          | Splitting scheme for the integration of Dissipative Particle Dynamics (DPD)                                                         | [  &dagger; &Dagger;](#simulation-mode) | [Doc](#dpd)      |
| `dpde_ser`               | Splitting with Energy Reinjection (SER) scheme for the integration of Dissipative Particle Dynamics with energy conservation (DPDE) | [  &dagger; &Dagger;](#simulation-mode) | [Doc](#dpde)     |

In order to increase the time and length scales that can be simulated, groups of atoms may be represented by a single particle. This effectively reduces the number of degrees of freedom which are explictly handled. The coarse-grained particles (or mesoparticles) are described by their positions \f$ q \f$, momenta \f$ p \f$ and possibly their internal energies \f$ \varepsilon_i \f$ taht takes into account all the coarse-grained degrees of freedom. The mesoparticles interact with an interaction [potential](#potential) \f$ U \f$ and with fluctuation/disipation forces. These methods require the keyword [`mode`](#simulation-mode) to be set to `single_spec_meso` or `single_spec_smooth`.

### <a name="dpd"></a> Dissipative Particle Dynamics

In DPD @cite hoogerbrugge_1992 @cite espanol_1995, mesoparticles interact with a potential and pairwise dissipation and stochastic forces which fix the temperature \f$ T \f$ in average. For a pair \f$ (i,j) \f$ of particles, we introduce their relative velocity \f$ v_{ij} = \frac{p_i}{m_i}-\frac{p_j}{m_j} \f$ and distance \f$ \require{AMSmath}r _{ij} = \left\lvert q_i-q_j \right\rvert \f$. The range of the fluctuation/dissipation interaction is limited by a cut-off function 
\f[
\chi(r) = \left(-\frac{1}{r_{\rm cut}}\right)^21_{r<r_{\rm cut}}.
\f]
with the same cut-off radius $r_{\rm cut}$ as the potential. The equations of motion for DPD read
\f{equation}{
\left\{
\begin{aligned}
dq_i & = \frac{p_i}{m_i}dt\\
dp_i &= -\nabla_{q_i} U(q)dt + \sum_{j\neq i} -\Gamma_{ij} \chi(r_{ij})v_{ij}dt + \sqrt{\chi(r_{ij})} \Sigma_{ij} dB_{ij},
\end{aligned}
\right.
\f}
where \f$ B_{ij} \f$ are 3-dimensional vectors of Brownian motions such that \f$ B_{ji} = - B_{ij} \f$.
The friction and fluctuation coefficients (respectively \f$ \Gamma_{ij} \f$ and \f$ \Sigma_{ij} \f$) are \f$ 3\times3 \f$ matrices under the form 
\f{equation}{
\begin{aligned}
\Gamma_{ij} $= \gamma^{\parallel}P_{ij}^{\parallel} + \gamma^{\perp}P_{ij}^{\perp},\\
\Sigma_{ij} $= \sigma^{\parallel}P_{ij}^{\parallel} + \sigma^{\perp}P_{ij}^{\perp},
\end{aligned}
\f}
where \f$ P_{ij}^{\parallel} \f$ is the projection matrix on the line of the centers of mass and \f$ P_{ij}^{\perp} \f$ in the orthogonal plane. The coefficients \f$ \gamma^{\parallel} \f$ and \f$ \gamma^{\perp} \f$ are fixed with the `meso_gamma_para` and `meso_gamma_perp` keywords (see [Mesoparticles](#mesoaprticle)). For \f$ \theta \in \{\parallel,\perp\} \f$, the friction and fluctuation coefficients are related with
\[
(\sigma^{\theta})^2 = 2\gamma^{\theta}k_{\rm B}T.
\]

The integration of DPD is performed thanks to a splitting strategy. The Hamiltonian part is handled with a Velocity-Verlet algorithm while the fluctuation/dissipation part is integrated with a Euler-Maruyama scheme. Starting from an initial configuration \f$ (q^n,p^n) \f$ and using a time step \f$ \Delta t \f$, the updated positions and  momenta at step \f$ n+1 \f$ are given by:
\f{equation}{
\left\{
\begin{aligned}
p_i^{n+\frac{1}{2}} &= p_i^n -\nabla_{q_i} U(q^n) \frac{\Delta t}{2},\\
q_i^{n+1} &= q_i^n + \frac{p_i^{n+\frac{1}{2}}}{m_i} \Delta t, \\
\tilde{p}_i^{n+1} &= p_i^{n+\frac{1}{2}} -\nabla_{q_i} U(q^{n+1}) \frac{\Delta t}{2},\\
p_i^{n+1} &= \tilde{p}_i^{n+1} + \sum_{j\neq i} -\Gamma_{ij}^{n+1} \chi(r_{ij}^{n+1})\tilde{v}_{ij}^{n+1} \Delta t + \sqrt{\Delta t\chi(r_{ij}^{n+1})} \Sigma_{ij}^{n+1} G_{ij},
\end{aligned}
\right.
\f}
with \f$ G_{ij} \f$ standard Gaussian variables.

| Required parameters      | Description                       | Unit         |
| ------------------------ | --------------------------------- | ------------ |
| `delta`                  | Time step \f$ \Delta t\f$         | s            |
| `temperature_thermostat` | Target temperature \f$ T \f$      | K            |
| `number_of_steps`        | Number of steps in the simulation | integer      |

### <a name="dpde"></a> Dissipative Particle Dynamics with Energy conservation

DPDE @cite avalos_1997 @cite espanol_1997 has been introduced to extend DPD to non isothermal situations. The coarse-grained internal degrees of freedom are now represented with an internal energy \f$ \varepsilon_i \f$ associated to each mesoparticle. An [equation of state](#eos) allows to compute an internal temperature \f$ T_i \f$ from the internal energy. The total energy of the system, which reads
\f[
E(q,p,\varepsilon) = U(q) + \sum_{i=1}^N \varepsilon_i + \frac{p_i^2}{2m},
\f]
is preserved by the dynamics. Compared to DPD, a dependence on the internal energies is introduced in the friction coefficients, now denoted by \f$ \gamma_{ij}^{\theta} \f$ for \f$ \theta \in \{\parallel,\perp\} \f$. 
The values of the fluctuation amplitude are determined from the `meso_gamma_para` and `meso_gamma_perp` keywords (denoted as \f$ \gamma^{\theta} \f$) with the relation
\f[
\sigma^{\theta} = 2\sqrt{\gamma^{\theta}k_{\rm B}T_{\rm ref}},
\f]
and \f$ T_{\rm ref} = 1 \f$ K.
The friction coefficients are given by
\f[
\gamma_{ij}^{\theta} = \frac{\sigma^{\theta}}{4}\left(\frac{1}{k_{\rm B}T_i} + \frac{1}{k_{\rm B}T_j}\right).
\f]
The variation in the kinetic energy generated by the fluctuation/dissipation forces is redistributed in the internal energies. Such a mechanism allows the internal and external degrees of freedom to exchange energy and ensures the equilibration of the internal and kinetic temperatures.
The equations of motion for DPDE read
\f{equation}{
\left\{
\begin{aligned}
dq_i & = \frac{p_i}{m_i}dt\\
dp_i &= -\nabla_{q_i} U(q)dt + \sum_{j\neq i} -\Gamma_{ij} \chi(r_{ij})v_{ij}dt + \sqrt{\chi(r_{ij})} \Sigma_{ij} dB_{ij},\\
d\varepsilon_i &= \frac{1}{2}\left[v_{ij}^T\Gamma_{ij}v_{ij} - \frac{1}{\mu_{ij}}{\rm Tr}(\Sigma_{ij}\Sigma_{ij}^T)\right]dt -\frac{1}{2} v_{ij}^T\Sigma_{ij}dB_{ij}.
\end{aligned}
\right.
\f}
where \f$ B_{ij} \f$ are 3-dimensional vectors of Brownian motions such that \f$ B_{ji} = - B_{ij} \f$ and \f$ \mu_{ij} = 2\left(\frac{1}{m_i}+\frac{1}{m_j}\right)^{-1} \f$ is the reduced mass.

The equations of motion of DPDE can be integrated with a Splitting with Energy Reinjection (SER) scheme @cite homman_2016. It consists in using a Velocity-Verlet scheme to discretize the Hamiltonian part of the dynamic. The fluctuation/dissipation part is then handled as follows. The momentum equation is integrated with an Euler-Maruyama scheme. Then the energy variation is redistributed in a symmetric pairwise form in the internal energies. This ensures the consevation of the energy. Starting from an initial configuration \f$ (q^n,p^n,\varepsilon^n) \f$ and using a time step \f$ \Delta t \f$, the SER scheme gives the positions, momenta and internal energies at step \f$ n+1 \f$ as:
\f{equation}{
\left\{
\begin{aligned}
p_i^{n+\frac{1}{2}} &= p_i^n -\nabla_{q_i} U(q^n) \frac{\Delta t}{2},\\
q_i^{n+1} &= q_i^n + \frac{p_i^{n+\frac{1}{2}}}{m_i} \Delta t, \\
\tilde{p}_i^{n+1} &= p_i^{n+\frac{1}{2}} -\nabla_{q_i} U(q^{n+1}) \frac{\Delta t}{2},\\
p_i^{n+1} &= \widetilde{p}_i^{n+1} + \sum_{j\neq i}\widetilde{\Gamma}_{ij}^{n+1}\widetilde{v}_{ij}^{n+1}\Delta t +\widetilde{\Sigma}_{ij}^{n+1}G_{ij}^n\sqrt{\Delta t},\\
\varepsilon_i^{n+1} &= \widetilde{\varepsilon}_i^{n+1} + \sum_{j\neq i}\frac{1}{2}(\widetilde{\Gamma}_{ij}^{n+1}\widetilde{v}_{ij}^{n+1}\Delta t - \widetilde{\Sigma}_{ij}^{n+1}G_{ij}^n\sqrt{\Delta t})\cdot\left(\widetilde{v}_{ij}^n+\frac{\delta p_i^n}{2m_i}-\frac{\delta p_j^n}{2m_j}\right),
\end{aligned}
\right.
\f}
where \f$ G_{ij}^n \f$ are independent standard Gaussian variables. The variation of momentum \f$ \delta p_i^n \f$ due to the fluctuation/dissipation part is given by
\f[
\delta p_i^n = \sum_{j\neq i} -\widetilde{\Gamma}_{ij}^{n+1}\widetilde{v}_{ij}^n\Delta t + \widetilde{\Sigma}_{ij}^nG_{ij}^n\sqrt{\Delta t}.
\f]

| Required parameters      | Description                       | Unit         |
| ------------------------ | --------------------------------- | ------------ |
| `delta`                  | Time step \f$ \Delta t\f$         | s            |
| `number_of_steps`        | Number of steps in the simulation | integer      |

## <a name="mesoparticle"></a> 3. Mesoparticle types

To define the types of particles involved into the simulation, several parameters need to be given. For each parameter, the associated input field reads an array allowing to specifiy to specify consecutively the data corresponding to several types of mesoparticles.

| Keyword            | Description                                                                                                 | Unit             | Modes                                   |
| ------------------ | ----------------------------------------------------------------------------------------------------------- | ---------------- | --------------------------------------- |
| `meso_names`       | Names of the mesoparticles                                                                                  | string, no space | [&dagger; &Dagger;](#Simulation mode) |
| `meso_masses`      | Masses of the mesoparticles                                                                                 | kg               | [&dagger; &Dagger;](#Simulation mode) |
| `meso_sizes`       | Number of atoms within a mesoparticle                                                                       | integer          | [&dagger; &Dagger;](#Simulation mode) |
| `meso_gamma_para`  | Friction coefficients \f$ \gamma^{\parallel} \f$ in the direction parallel to the line of centers of mass   | s\f$^{-1}\f$     | [&dagger; &Dagger;](#Simulation mode) |
| `meso_gamma_ortho` | Friction coefficients \f$ \gamma^{\perp} \f$ in the direction orthogonal to the line of centers of mass     | s\f$^{-1}\f$     | [&dagger; &Dagger;](#Simulation mode) |

## <a name="interaction"></a> 4. Interactions

### <a name="potential"></a> Potential functions
In coarse-grained models, potential functions describe the interactions between particles, allowing to compute the potential energy and forces. For each type of potential (Lennard-Jones, EAM,...), each keyword is an array which are filled with the input data for each interaction.

#### <a name="lennard-jones"></a> Lennard-Jones

The Lennard-Jones potential involves two types of atoms (specified with `lj_types_A` and `lj_types_B`). For a pair of particles \f$ (i,j) \f$, with \f$i\f$ of type \f$A\f$ and \f$j\f$ of type \f$B\f$, the interaction energy depends only on the distance \f$ \require{AMSmath} r_{ij} = \left\lvert q_i-q_j \right\rvert \f$ and is given by
\f[
U_{\rm LJ}(r_{ij}) = 4\varepsilon\left(\left(\frac{\sigma}{r_{ij}}\right)^{12}-\left(\frac{\sigma}{r_{ij}}\right)^6\right),
\f]
with two parameters: \f$ \sigma \f$, such that \f$U_{\rm LJ}'(2^{\frac{1}{6}}\sigma)=0 \f$, controls the distance at which the potential energy is mjnimum while \f$ \varepsilon \f$ fixes the depth of the energy well.
The potential is actually set to \f$ 0 \f$ if \f$ r_{ij} \f$ is larger than a cut-off radius \f$ r_{\rm cut} \f$.

| Keyword      | Description                                  | Unit             | 
| ------------ | -------------------------------------------- | ---------------- | 
| `lj_type_A`  | Names of the particles                       | string, no space | 
| `lj_type_B`  | Names of the particles                       | string, no space | 
| `lj_rcut`    | Cut-off radius                               | m                | 
| `lj_epsilon` | Depth of the energy well \f$ \varepsilon \f$ | J                | 
| `lj_sigma`   | Interaction length \f$ \sigma \f$            | m                | 

### <a name="fd"></a> Fluctuation/dissipation

On top of the potential interactions, DPD and DPDE includes a fluctuation/dissipation part which is activated with the following keywords.

### <a name="dpd-interaction"></a> DPD 

| Keyword       | Description            | Unit             | 
| ------------- | ---------------------- | ---------------- | 
| `dpd_type_A`  | Names of the particles | string, no space | 
| `dpd_type_B`  | Names of the particles | string, no space | 

### <a name="dpde-interaction"></a> DPDE 

| Keyword        | Description            | Unit             | 
| -------------- | ---------------------- | ---------------- | 
| `dpde_type_A`  | Names of the particles | string, no space | 
| `dpde_type_B`  | Names of the particles | string, no space | 

## <a name="eos"></a> 5. Equations of state

When using DPDE, an equation of state allows to compute the internal temperature \f$ T_i \f$ of a mesoparticle as a function of the internal energy \f$ \varepsilon_i \f$. 

### <a name="eos-dpde"></a> Constant heat capacity (CHC)

The simplest model consists in assuming a constant heat capacity \f$ C_V \f$. The temperature is then given by
\f[
T_i = \frac{\varepsilon_i}{C_V}.
\f]

| Keyword     | Description                                          | Unit             |
| ----------- | ---------------------------------------------------- | ---------------- |
| `dpde_type` | Name of the particles with the CHC equation of state | string, no space |
| `dpde_cv`   | Heat capcity \f$ C_v \f$                             | J.K\f$^{-1}\f$   |


### <a name="eos-sdpd"></a> Other equations of state

More complex equations of state are implemented in ExaStamp to be used with SDPD and relate the entropy with the energy and the density. If they are used in the context of DPDE, they resort to a local estimation of the density in its Verlet cell (computed as the total mass within the cell divided by its volume). We refer to the [SDPD help page](md_markdown_sdpd.html#eos) for their documentation.

## <a name="domain"></a> 6. Setting the domain

In order to initialize the simulation box, the following parameters are required. Note that the `extension` keyword is overridden if the particle configuration is initialized from a dump file or from a lattice configuration (see [Particle configuration](#configuration)).

| Keyword      | Description                                  | Unit             | 
| ------------ | -------------------------------------------- | ---------------- | 
| `origin`     | Origin of the domain                         | m (3-dim array)  | 
| `extension`  | Extension of the domain                      | m (3-dim array)  | 

### <a name="boundary"></a> Boundary conditions

The boundary conditions are specified thanks to the `boundary_conditions` keyword. An array of three values chosen in the following table are expected (one for each direction). 

| Keyword: `boundary_conditions` | Description                                                                                   |
| ------------------------------ | --------------------------------------------------------------------------------------------- |
| `free`                         | No boundary conditions. Particles may leave the domain (see [Expanding the domain](#expand)). |
| `periodic`                     | Periodic boundary conditions.                                                                 |
| `wall`                         | A wall is located at each end of the domain (see [Walls](#wall)).                             |
| `free_wall`                    | Free boundary conditions at the lower end, wall at the upper end.                             |
| `wall_free`                    | Wall at the lower end, free boundary conditions at the upper end.                             |

### <a name="expand"> Expanding the domain

If free boundary conditions are used, it is possible to choose whether the particles leaving the domain should be lost or whether the simulation box is extented to keep track of them. If the keyword `expand_free_limits` is set to true, the simulation is stopped when a particle is about to leave the domain. A dump file is then generated and the extension of the domain is multiplied by a factor set with the keyword `expand_factor` in the problematic direction.

| Keyword              | Description                                                                      | Unit    | 
| -------------------- | -------------------------------------------------------------------------------- | ------- | 
| `expand_free_limits` | If set to true, activates the expansion of the domain when a particle leaves it. | boolean | 
| `expand_factor`      | Factor by which the domain is extended                                           | double  | 

### <a name="configuration"> Initial particle configuration

The initial particle configuration can either be generated on a given lattice or imported from a dump file. The keyword `initialization` is responsible for this choice.

| Keyword: `initialization` | Description                                                                  |
| ------------------------- | ---------------------------------------------------------------------------- |
| `default`                 | Initialization from a [lattice](#lattice)                                    |
| `legacy`                  | Initialization from a [MpiIO dump file](#mpiio-read) (compatible with Stamp) |
| `hercule`                 | Initialization from a [Hercule database](#hercule-read)                      |

#### <a name="lattice"></a> Initializing a lattice

In order to create an initial lattice, the type of the lattice must be chose with the keyword `lattice_type`. This fixes \f$ n_{\rm at} \f$ which is the number of atoms in a lattice cell.

| Keyword: `lattice_type` | Description                                                |
| ----------------------- | ---------------------------------------------------------- |
| `bcc`                   | Body centered cubic lattice (\f$ n_{\rm at} = 2 \f$ atoms) |
| `fcc`                   | Face centered cubic lattice (\f$ n_{\rm at} = 4 \f$ atoms) |
| `sc`                    | Simple cubic lattice (\f$ n_{\rm at} = 1\f$ atom)          |
| `diam100`               | Diamond lattice (\f$ n_{\rm at} = 8\f$ atoms)              |

The positions of the atoms are then initialized on the sites of the lattice cell.

| Keyword             | Description                                                                                   | Unit                         | 
| ------------------- | --------------------------------------------------------------------------------------------- | ---------------------------- | 
| `lattice_parameter` | Dimension of a lattice cell                                                                   | m                            |
| `lattice_atoms`     | Array of size \f$ n_{\rm at} \f$ with the atom types associated with each site of the lattice | string (x\f$ n_{\rm at} \f$) |
| `n_cells`           | Number of repetition of the lattice cell in each direction                                    | integer (3-dim array)        |

The initial velocities are distributed along a Maxwellian with the initial tempertaure set thanks to the `init_temperature` keyword. The initial internal energies are chosen such that the internal temperature for each particle is equal to the value of `init_tint`.

| Keyword            | Description                  | Unit | 
| ------------------ | ---------------------------- | ---- | 
| `init_temperature` | Initial temperature          | K    |
| `init_tint`        | Initial internal temperature | K    |

#### <a name="mpiio-read"></a> Reading a MpiIO dump file

One of the two following keywords need to be set (`particle_file` has the priority).

| Keyword          | Description                                                                    | Unit    | 
| ---------------- | ------------------------------------------------------------------------------ | ------- | 
| `init_step`      | Initializes from `StampV3prot_xxxxxxxxx.MpiIO` corresponding to the given step | integer |
| `particles_file` | Initializes from the given MpiIO file                                          | string  |

#### <a name="hercule-read"></a> Reading a Hercule database

| Keyword          | Description                                                    | Unit    | 
| ---------------- | -------------------------------------------------------------- | ------- | 
| `init_step`      | Initializes from a Hercule database starting at the given step | integer |

## <a name="output"></a> 7. Output

ExaStamp provide different types of output: global data (*e.g.* total energy, temperature,...) printed in a log file, particle data (*e.g.* positions,...), cell data (local averages like temperature or pressure). In order to restart simulations, dump file may also be generated. The table below summarizes the output options.

| Keyword        | Description                                                                                  |
| -------------- | -------------------------------------------------------------------------------------------- |
| `log_rate`     | Number of steps between each [log](#log) output (-1 to disable)                              |
| `output_rate`  | Number of steps between each particle output (-1 to disable)                                 |
| `output_type`  | Type of particle output ([dat](#output-dat), [vtk](#output-vtk), [hercule](#output-hercule)) |
| `output_cells` | Activate a [vtk](#output-vtk) output of cell data                                            |
| `output_dir`   | Directory where the particle output is stored (default is current directory)                 |
| `dump_rate`    | Number of steps between each dump (-1 to disable)                                            |
| `dump_type`    | Type of dump ([legacy](#dump-legacy), [hercule](#dump-hercule))                              |
| `dump_dir`     | Directory where the dump is stored (deault is current directory)                             |

### <a name="log"></a> Log file

The log file shows the iteration number, the number of particles, the total energy, the potential energy, the internal energy, the chemical energy, the temperature, the pressure and the computation time per particle per task per iteration.

### <a name="output-dat"></a> Particle output with DAT files

If the keyword `output_type` is set to `dat`, the particle output will be generated under the form of a single text file (.DAT) for each timestep. It contains the positions, ids, types, velocities, internal energies and progress variables of each particle. 

### <a name="output-hercule"></a> Particle and cell output with VTK files

If the keyword `output_type` is set to `vtk`, the particle output will be generated under the form of a single VTK file for each timestep. It contains the positions, ids, types, atomic numbers, velocities, internal energies and progress variables of each particle. It may be read with the ParaView software with the PointSprite plugin enabled.

If the keyword `output_cells` is set to `true`, a cell output will be generated under the form of a single VTK file for each timestep, independently of the value of `output_type`. It contains the positions, pressure tensors, temperatures and densities of each Verlet cell. It may be read with the ParaView software.

### <a name="output-hercule"></a> Particle and cell output with a Hercule database

If the keyword `output_type` is set to `hercule`, the particle output will be generated under the form of a single Hercule database for all timesteps (which may be split in several files if they become too large). It contains the positions, ids, types, atomic numbers, velocities, internal energies and progress variables of each particle. It may be read with the Love software (CEA version of ParaView).

### <a name="dump-legacy"></a> Dump with MpiIO files

If the keyword `dump_type` is set to `legacy`, the dump file will be generated under the form of a single MpiIO file for each timestep. It should be compatible with Stamp.

### <a name="dump-hercule"></a> Dump with a Hercule database

If the keyword `dump_type` is set to `hercule`, the dump file will be generated under the form of a Hercule database for all timesteps (which may be split in several files if they become too large).

## <a name="parallel"></a> 8. Parallelisation

| Keyword                | Description                               |
| ---------------------- | ----------------------------------------- |
| `decoupage`            | Number of MPI processes in each direction |
| `max_threads_per_node` | Maximum number of TBB threads per process |

# <a name="example"></a> Examples

## Input file for DPD
@include ../../tests/single-mat-dpd.xsp

## Input file for DPDE
@include ../../tests/single-mat-dpde.xsp


