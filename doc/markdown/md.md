# Molecular Dynamics

The main purpose of ExaStamp is to allow for the efficient simulation of Molecular Dynamics (MD). This help page summarizes the keywords needed to set up such a simulation.

## <a name="simulation-mode"></a> 1. Simulation mode 

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
| `newton_verlet_leapfrog` | Leapfrog Verlet scheme for the integration of Hamiltonian dynamics (NVE)                                                            | [* &dagger; &Dagger;](#Simulation mode) | [Doc](#leapfrog) |
| `newton_verlet_velocity` | Velocity Verlet scheme for the integration of Hamiltonian dynamics (NVE)                                                            | [* &dagger; &Dagger;](#simulation-mode) | [Doc](#velocity) |
| `langevin_splitting`     | Integration of the Langevin dynamics (NVT)                                                                                          | [* &dagger; &Dagger;](#simulation-mode) | [Doc](#langevin) |

In classical MD, the atoms are described by their position \f$ q \f$ and velocities \f$ p \f$. The interactions between atoms are represented with a [potential](#potential) function \f$ U \f$ that depends on the atomic positions. For a system of \f$ N \f$ atoms indexed by \f$ i \f$ and of mass \f$ m_i \f$, the total energy of the system is defined as the sum of the kinetic and potential energies as
\f[
H(q,p) = \sum_{i=1}^N \frac{p_i^2}{2m_i} + U(q).
\f]
Any [`mode`](#simulation-mode) allows to run MD simulations in ExaStamp.

### <a name="nve"></a> NVE ensemble

In the NVE ensemble , the number of particles \f$ N \f$, the volume \f$ V \f$, and the total energy \f$ E \f$ are kept constant throughout the simulation. This is achieved by integrating a Hamiltonian dynamics:
\f{equation}{
\left\{
\begin{aligned}	
dq_i &= \frac{p_i}{m_i}dt, \\
dp_i &= -\nabla_{q_i} U(q)dt.
\end{aligned}
\right.
\f}

To ensure a long-term conservation of the energy \f$ H(q,p) \f$, it is required to integrate the equations of motion with a symplectic scheme. ExaStamp implememts two such schemes that are both versions of the Verlet algorithm. The step of the time discretization is denoted by \f$ \Delta t \f$.

#### <a name="leapfrog"></a> Leapfrog Verlet

Starting from an initial configuration \f$ (q^n,p^n) \f$, the Leapfrog Verlet scheme gives the positions and momenta at step \f$ n+1 \f$ as:
\f{equation}{
\left\{
\begin{aligned}
q_i^{n+\frac{1}{2}} &= q_i^n + \frac{p_i^n}{m_i} \frac{\Delta t}{2}, \\
p_i^{n+1} &= p_i^n -\nabla_{q_i} U(q^{n+\frac{1}{2}}) \Delta t,\\
q_i^{n+1} &= q_i^{n+\frac{1}{2}} + \frac{p_i^{n+1}}{m_i} \frac{\Delta t}{2}.
\end{aligned}
\right.
\f}
The scheme is activated by setting the keyword [`scheme`](#scheme) to `newton_verlet_leapfrog`.

| Required parameters | Description                       | Unit    |
| ------------------- | --------------------------------- | ------- |
| `delta`             | Time step \f$ \Delta t \f$        | s       |
| `number_of_steps`   | Number of steps in the simulation | integer |

#### <a name="velocity"></a> Velocity Verlet

Starting from an initial configuration \f$ (q^n,p^n) \f$, the Velocity Verlet scheme gives the positions and momenta at step \f$ n+1 \f$ as:
\f{equation}{
\left\{
\begin{aligned}
p_i^{n+\frac{1}{2}} &= p_i^n -\nabla_{q_i} U(q^n) \frac{\Delta t}{2},\\
q_i^{n+1} &= q_i^n + \frac{p_i^{n+\frac{1}{2}}}{m_i} \Delta t, \\
p_i^{n+1} &= p_i^{n+\frac{1}{2}} -\nabla_{q_i} U(q^{n+1}) \frac{\Delta t}{2}.\\
\end{aligned}
\right.
\f}
The scheme is activated by setting the keyword [`scheme`](#scheme) to `newton_verlet_velocity`.

| Required parameters | Description                       | Unit    |
| ------------------- | --------------------------------- | ------- |
| `delta`             | Time step \f$ \Delta t \f$        | s       |
| `number_of_steps`   | Number of steps in the simulation | integer |

### <a name="nvt"></a> NVT ensemble

In the NVT ensemble , the number of particles \f$ N \f$, the volume \f$ V \f$, and the temperature \f$ T \f$ are preserved throughout the simulation. A thermostat is used to keep the temperature at the desired value.

#### <a name="langevin"></a> Langevin thermostat

The Langevin dynamics achieves a constant temperature simulation (in average) through the use of a stochastic thermostat. An additional dissipation force along with a random fluctation are added to the Hamiltonian dynamics and related together by a fluctuation-dissipation relation. With \f$ T \f$ the desired temperature, \f$ \gamma \f$ the friction parameter and \f$ B_t \f$ a Brownian motion, the equations of motion for the Langevin dynamics read
\f{equation}{
\left\{
\begin{aligned}	
dq_i &= \frac{p_i}{m_i}dt, \\
dp_i &= -\nabla_{q_i} U(q)dt -\gamma p_idt + \sqrt{2m_i\gamma k_{\rm B}T} dB_t.
\end{aligned}
\right.
\f}
Here \f$ k_{\rm B} \f$ is the Boltzmann constant.

The Langevin dynamics is discretized with a splitting strategy. The Hamiltonian dynamics is integrated with the [Velocity Verlet](#langevin) scheme while the fluctuation and dissipation terms are integrated analytically as an Ornstein-Uhlenbeck process. To that end, we introduce standard Gaussian variables \f$ G_i \f$ along with the quantities
\f{equation}{
\left\{
\begin{aligned}
\alpha_{\Delta t} &= \exp\left(-\gamma\Delta t\right),\\
\zeta_{\Delta t} &= \sqrt{k_{\rm B}T\left(1-(\alpha_{\Delta t})^2\right)}.
\end{aligned}
\right.
\f}
The overall scheme starting from an initial configuration \f$ (q^n,p^n) \f$ reads
\f{equation}{
\left\{
\begin{aligned}
p_i^{n+\frac{1}{2}} &= p_i^n -\nabla_{q_i} U(q^n) \frac{\Delta t}{2},\\
q_i^{n+1} &= q_i^n + \frac{p_i^{n+\frac{1}{2}}}{m_i} \Delta t, \\
\tilde{p}_i^{n+1} &= p_i^{n+\frac{1}{2}} -\nabla_{q_i} U(q^{n+1}) \frac{\Delta t}{2},\\
p_i^{n+1} &= \alpha_{\Delta t}\tilde{p}_i^{n+1} + \sqrt{m_i}\zeta_{\Delta t}G_i^n.
\end{aligned}
\right.
\f}
The scheme is activated by setting the keyword [`scheme`](#scheme) to `langevin_splitting`.

| Required parameters      | Description                       | Unit         |
| ------------------------ | --------------------------------- | ------------ |
| `delta`                  | Time step \f$ \Delta t\f$         | s            |
| `friction_langevin`      | Friction parameter \f$ \gamma \f$ | s\f$^{-1}\f$ |
| `temperature_thermostat` | Target temperature \f$ T \f$      | K            |
| `number_of_steps`        | Number of steps in the simulation | integer      |

## <a name="atom"></a> 3. Atom types

To define the types of atom involved into the simulation, several parameters need to be given. For each parameter, the associated input field reads an array allowing to specifiy to specify consecutively the data corresponding to several types of atom.

| Keyword        | Description          | Unit             | Modes                                   |
| -------------- | -------------------- | ---------------- | --------------------------------------- |
| `atom_names`   | Names of the atoms   | string, no space | [* &dagger; &Dagger;](#Simulation mode) |
| `atom_z`       | Atomic numbers       | integer          | [* &dagger; &Dagger;](#Simulation mode) |
| `atom_masses`  | Masses of the atoms  | kg               | [* &dagger; &Dagger;](#Simulation mode) |
| `atom_charges` | Charges of the atoms | C                | [* &dagger; &Dagger;](#Simulation mode) |

An example in a simulation involving only copper and zinc reads:

```
atom_names   = [      Copper,        Zinc]
atom_z       = [          29,          30]
atom_masses  = [ 105.486E-27, 108.531E-27]
atom_charges = [           0,           0]
```

## <a name="potential"></a> 4. Potential functions

In classical MD, potential functions describe the interactions between particles, allowing to compute the potential energy and forces. For each type of potential (Lennard-Jones, EAM,...), each keyword is an array which are filled with the input data for each interaction. For instance, in a system with copper and zinc interacting with Lennard-Jones potentials, the input file may look like:

```
lj_type_A  = [   Copper,      Zinc,    Copper]
lj_type_B  = [   Copper,      Zinc,      Zinc]
lj_rcut    = [0.568E-09, 0.610E-09, 0.589E-09]
lj_epsilon = [9.340E-20, 2.522E-20, 4.853E-20]
lj_sigma   = [0.227E-09, 0.244E-09, 0.236E-09]
```

### <a name="lennard-jones"></a> Lennard-Jones

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

### <a name="sutton-chen"></a> Sutton-Chen

### <a name="vniitf"></a> VNIITF

### <a name="meam"></a> MEAM 

## <a name="domain"></a> 5. Setting the domain

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

The initial velocities are distributed along a Maxwellian with the initial tempertaure set thanks to the `init_temperature` keyword.

| Keyword            | Description         | Unit | 
| ------------------ | ------------------- | ---- | 
| `init_temperature` | Initial temperature | K    |

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

## <a name="output"></a> 6. Output
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

## <a name="parallel"></a> 7. Parallelisation

| Keyword                | Description                               |
| ---------------------- | ----------------------------------------- |
| `decoupage`            | Number of MPI processes in each direction |
| `max_threads_per_node` | Maximum number of TBB threads per process |


# <a name="example"></a> Examples

## Input file for Lennard-Jones (NVE)
@include ../../tests/multi-mat-lj.xsp

## Input file for Sutton-Chen (NVE)
@include ../../tests/single-mat-sc.xsp

## Input file for VNIITF (NVE)
@include ../../tests/single-mat-vniitf.xsp

## Input file for MEAM (NVE)
@include ../../tests/single-mat-meam.xsp
