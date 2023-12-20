# Smoothed Dissipative Particle Dynamics

Smoothed Dissipative Particle Dynamics is a mesoscopic model coupling a particle discretization of the Navier-Stokes eqautions and thermal fluctuations of microscopic origin. It has been included in ExaStamp.

## <a name="simulation-mode"></a> 1. Simulation mode 

The choice of a simulation mode is mandatory in ExaStamp. Each enables a different set of functionalities. The mode is selected by adding the `mode` keyword with the desired value. In the rest of the documentation, the options requiring a specific mode are indicated and signaled with the symbols associated with the compatible modes.

| Keyword: `mode`      | Description                                            | Symbol   |
| -------------------- | ------------------------------------------------------ | :------: |
| `single_spec_atom`   | Used for classical MD                                  | *        |
| `single_spec_meso`   | Used for coarse-grained dynamics (DPD, DPDE)           | &dagger; |
| `single_spec_smooth` | Used for Smoothed Dissipative Particle Dynamics (SDPD) | &Dagger; |


## <a name="scheme"></a> 2. Scheme 

The choice of a numerical scheme is mandatory. It tells ExaStamp what dynamics to run. The modes compatible with each scheme are indicated in the following table and links to a more thorough documentation on the dynamics and their integration scheme is given. The schemes available to perform classical MD simulations are listed in the following table.

| Keyword: `scheme`   | Description                                                                                                                         | Modes                                   | Link             |
| ------------------- | -----------------------------------------------------------------------------------------------------------------------------------------------------  | :-------------------------------------: | :--------------: | 
| `sdpd_vv_ser`       | Splitting with Energy Reinjection (SER) scheme for the integration of Smoothed Dissipative Particle Dynamics (SDPD)                                   | [&Dagger;](#simulation-mode) | [Doc](#sdpd)     |
| `sdpd_langevin_ser` | Splitting with Energy Reinjection (SER) scheme for the integration of Smoothed Dissipative Particle Dynamics (SDPD) coupled with a Langevin thermostat | [&Dagger;](#simulation-mode) | [Doc](#sdpd)     |

### <a name="sdpd"></a> Smoothed Dissipative Particle Dynamics

SDPD @cite espanol_2003 couples a particle discretization of the Navier-Stokes equations (see SPH @cite lucy_1977 @cite  monaghan_1977) with thermal fluctuations varying consistently with the resolution. The domain is discretized in fluid particles (or smooth particles) representing a mass \f$ m_i \f$ of the fluid and acting as interpolation nodes for the resolution of the hydrodynamics equations through a [smoothing kernel](#kernel) \f$ W \f$. The smooth particles are described by their position \f$ q_i \f$, momenta \f$ p_i \f$ and internal energies \f$ \varepsilon_i \f$. A density \f$ \rho_i \f$ is associated to each particle as an average over neighboring particles:
\f{equation}{
\rho_i = \sum_{i=1}^N m_j W(r_{ij}),
\f}
depending on the inter-particle distance \f$ \require{AMSmath} r_{ij} = \left\lvert q_i-q_j \right\rvert \f$.
The thermodynamic state (entropy \f$ S_i \f$, temperature \f$ T_i \f$, pressure \f$ P_i \f$ and heat capacity \f$ C_i \f$) of a fluid particle is determined from the density and the internal energy thanks to an [equation of state](#eos).
For a pair of particle \f$ (i,j) \f$, we also introduce \f$ F_{ij} \f$ defined such that
\f[
\nabla_{q_i}W(r_{ij}) = -F_{ij}\times(q_i-q_j),
\f]
along with their reduced mass \f$ \mu_{ij} = 2\left(\frac{1}{m_i}+\frac{1}{m_j}\right)^{-1} \f$ and their relative velocity \f$ v_{ij} = \frac{p_i}{m_i}-\frac{p_j}{m_j} \f$.
The total energy in a SDPD system is the sum of a kinetic term and of an internal term:
\f[
E(q,p,\varepsilon) = \sum_{i=1}^N \varepsilon_i + \frac{p_i^2}{2m}.
\f]

The equations of motion in the position \f$ q_i \f$, momentum \f$ p_i \f$ and internal energy \f$ \varepsilon_i \f$ variables, as reformulated in @cite faure_2016, read
\f{equation}{
\left\{
\begin{aligned}
dq_i &= \frac{p_i}{m_i}dt,\\
dp_i &= \sum_{j\neq i} m_im_j\left(\frac{P_i}{\rho_i^2}+\frac{P_j}{\rho_j^2}\right)F_{ij}(q_i-q_j)dt - \Gamma_{ij}v_{ij}dt + \Sigma_{ij}dB_{ij},\\
d\varepsilon_i &= \sum_{j\neq i} -m_im_j\frac{P_i}{\rho_i^2}F_{ij}(q_i-q_j)^Tv_{ij}dt +\frac{1}{2}\left[v_{ij}^T\Gamma_{ij}v_{ij} - \frac{1}{\mu_{ij}}{\rm Tr}(\Sigma_{ij}\Sigma_{ij}^T)\right]dt -\frac{1}{2} v_{ij}^T\Sigma_{ij}dB_{ij}.
\end{aligned}
\right.
\f}
They involve 3-dimensional vectors of Brownian motions \f$ B_{ij} \f$ such that \f$ B_{ji} = - B_{ij} \f$.
The friction and fluctuation coefficients (respectively \f$ \Gamma_{ij} \f$ and \f$ \Sigma_{ij} \f$) are \f$ 3\times3 \f$ matrices under the form 
\f{equation}{
\begin{aligned}
\Gamma_{ij} $= \gamma_{ij}^{\parallel}P_{ij}^{\parallel} + \gamma_{ij}^{\perp}P_{ij}^{\perp},\\
\Sigma_{ij} $= \sigma_{ij}^{\parallel}P_{ij}^{\parallel} + \sigma_{ij}^{\perp}P_{ij}^{\perp},
\end{aligned}
\f}
where \f$ P_{ij}^{\parallel} \f$ is the projection matrix on the line of the centers of mass and \f$ P_{ij}^{\perp} \f$ in the orthogonal plane. The coefficients \f$ \gamma_{ij}^{\theta} \f$ and \f$ \sigma_{ij}^{\theta} \f$ are related for \f$ \theta \in \{\parallel,\perp\} \f$ and depend on the bulk viscosity \f$ \eta \f$ and the shear viscosity \f$ \zeta \f$.
Introducing a cut-off function
\f[
\chi_{ij} = \frac{m_im_j}{\rho_i\rho_j}F_{ij},
\f]
the friction and fluctuation coefficients are given by
\f{equation}{
\begin{aligned}
d_{ij} &= k_{\rm B}\frac{T_iT_j}{(T_i+T_j)^2}\left(\frac{1}{C_i}+\frac{1}{C_j}\right),\\
\gamma_{ij}^{\parallel} &= \left(\frac{10}{3}\eta + 4\zeta\right)\chi_{ij}(1-d_{ij}),\\
\gamma_{ij}^{\perp} &= \left(\frac{5}{3}\eta-\zeta\right)\chi_{ij}(1-d_{ij}),\\
\sigma_{ij}^{\theta} &= 2\sqrt{\frac{\gamma^{\theta}}{1-d_ij}k_{\rm B}\frac{T_iT_j}{T_i+T_j}}.
\end{aligned}
\f}

#### <a name="sdpd-ser"> Splitting with Energy Reinjection

The equations of motion of SDPD can be integrated with a Splitting with Energy Reinjection (SER) scheme @cite homman_2016. It consists in using a Velocity-Verlet scheme to discretize the conservative part of the dynamic. The fluctuation/dissipation part is then handled as follows. The momentum equation is integrated with an Euler-Maruyama scheme. Then the energy variation is redistributed in a symmetric pairwise form in the internal energies. This ensures the consevation of the energy. Starting from an initial configuration \f$ (q^n,p^n,\varepsilon^n) \f$ and using a time step \f$ \Delta t \f$, the SER scheme gives the positions, momenta and internal energies at step \f$ n+1 \f$ as:
\f{equation}{
\left\{
\begin{aligned}
p_i^{n+\frac{1}{2}} &= p_i^n + \sum_{j\neq i} m_im_j\left(\frac{P_i^n}{(\rho_i^n)^2}+\frac{P_j^n}{(\rho_j^n)^2}\right)F_{ij}^n(q_i^n-q_j^n)\Delta t,\\
q_i^{n+1} &= q_i^n  + \frac{p_i^{n+\frac{1}{2}}}{m_i} \Delta t, \\
\widetilde{\varepsilon}_i^{n+1} &= \mathcal{E}(S_i^n,\rho_i^{n+1}),\\
\widetilde{p}_i^{n+1} &= p_i^{n+\frac{1}{2}} + \sum_{j\neq i} m_im_j\left(\frac{\widetilde{P}_i^{n+1}}{(\rho_i^{n+1})^2}+\frac{\widetilde{P}_j^{n+1}}{(\rho_j^{n+1})^2}\right)F_{ij}^{n+1}(q_i^{n+1}-q_j^{n+1})\Delta t,\\ 
p_i^{n+1} &= \widetilde{p}_i^{n+1} + \sum_{j\neq i}\widetilde{\Gamma}_{ij}^{n+1}\widetilde{v}_{ij}^{n+1}\Delta t +\widetilde{\Sigma}_{ij}^{n+1}G_{ij}^n\sqrt{\Delta t},\\
\varepsilon_i^{n+1} &= \widetilde{\varepsilon}_i^{n+1} + \sum_{j\neq i}\frac{1}{2}(\widetilde{\Gamma}_{ij}^{n+1}\widetilde{v}_{ij}^{n+1}\Delta t - \widetilde{\Sigma}_{ij}^{n+1}G_{ij}^n\sqrt{\Delta t})\cdot\left(\widetilde{v}_{ij}^n+\frac{\delta p_i^n}{2m_i}-\frac{\delta p_j^n}{2m_j}\right),
\end{aligned}
\right.
\f}
where \f$ G_{ij}^n \f$ are independent standard Gaussian variables and \f$ \mathcal{E} \f$ is the equation of state relating the internal energy \f$ \varepsilon_i \f$ with the entropy \f$ S_i \f$ and the density \f$ \rho_i \f$. The variation of momentum \f$ \delta p_i^n \f$ due to the fluctuation/dissipation part is given by
\f[
\delta p_i^n = \sum_{j\neq i} -\widetilde{\Gamma}_{ij}^{n+1}\widetilde{v}_{ij}^n\Delta t + \widetilde{\Sigma}_{ij}^nG_{ij}^n\sqrt{\Delta t}.
\f]


| Required parameters      | Description                       | Unit         |
| ------------------------ | --------------------------------- | ------------ |
| `delta`                  | Time step \f$ \Delta t\f$         | s            |
| `number_of_steps`        | Number of steps in the simulation | integer      |

#### <a name="sdpd-langevin"></a> SER with Langevin thermostat

This algorithm couples SDPD with a Langevin thermostat ensuring a constant temperature (in average). SDPD is integrated with the SER scheme while the Langevin dynamics is integrated analytically as an Ornstein-Uhlenbeck process.
Following the expressions of the [SER](#sdpd-ser) scheme and introducing independent standard Gaussian variables \f$ \overline{G}_i^n \f$, the updated configuration with the SEr-Langevin coupling reads:
\f{equation}{
\left\{
\begin{aligned}
p_i^{n+\frac{1}{2}} &= p_i^n + \sum_{j\neq i} m_im_j\left(\frac{P_i^n}{(\rho_i^n)^2}+\frac{P_j^n}{(\rho_j^n)^2}\right)F_{ij}^n(q_i^n-q_j^n)\Delta t,\\
q_i^{n+1} &= q_i^n  + \frac{p_i^{n+\frac{1}{2}}}{m_i} \Delta t, \\
\widetilde{\varepsilon}_i^{n+1} &= \mathcal{E}(S_i^n,\rho_i^{n+1}),\\
\overline{p}_i^{n+1} &= p_i^{n+\frac{1}{2}} + \sum_{j\neq i} m_im_j\left(\frac{\widetilde{P}_i^{n+1}}{(\rho_i^{n+1})^2}+\frac{\widetilde{P}_j^{n+1}}{(\rho_j^{n+1})^2}\right)F_{ij}^{n+1}(q_i^{n+1}-q_j^{n+1})\Delta t,\\ 
\widetilde{p}_i^{n+1} &= \alpha_{\Delta t}\overline{p}_i^{n+1} + \sqrt{m_i}\zeta_{\Delta t}\overline{G}_i^n,\\
p_i^{n+1} &= \widetilde{p}_i^{n+1} + \sum_{j\neq i}\widetilde{\Gamma}_{ij}^{n+1}\widetilde{v}_{ij}^{n+1}\Delta t +\widetilde{\Sigma}_{ij}^{n+1}G_{ij}^n\sqrt{\Delta t},\\
\varepsilon_i^{n+1} &= \widetilde{\varepsilon}_i^{n+1} + \sum_{j\neq i}\frac{1}{2}(\widetilde{\Gamma}_{ij}^{n+1}\widetilde{v}_{ij}^{n+1}\Delta t - \widetilde{\Sigma}_{ij}^{n+1}G_{ij}^n\sqrt{\Delta t})\cdot\left(\widetilde{v}_{ij}^n+\frac{\delta p_i^n}{2m_i}-\frac{\delta p_j^n}{2m_j}\right).
\end{aligned}
\right.
\f}
The integration of the Langevin dynamics involves the quantities:
\f{equation}{
\left\{
\begin{aligned}
\alpha_{\Delta t} &= \exp\left(-\overline{\gamma}\Delta t\right),\\
\zeta_{\Delta t} &= \sqrt{k_{\rm B}T\left(1-(\alpha_{\Delta t})^2\right)},
\end{aligned}
\right.
\f}
where \f$ T \f$ is the target temperature and \f$ \overline{\gamma} \f$ the friction parameter.

| Required parameters      | Description                                  | Unit         |
| ------------------------ | -------------------------------------------- | ------------ |
| `delta`                  | Time step \f$ \Delta t\f$                    | s            |
| `friction_langevin`      | Friction parameter \f$ \overline{\gamma} \f$ | s\f$^{-1}\f$ |
| `temperature_thermostat` | Target temperature \f$ T \f$                 | K            |
| `number_of_steps`        | Number of steps in the simulation            | integer      |

## <a name="smoothparticle"></a> 3. Smooth particles

To define the types of particles involved into the simulation, several parameters need to be given. For each parameter, the associated input field reads an array allowing to specifiy to specify consecutively the data corresponding to several types of smooth particles. The smoothing length \f$ h \f$ is determined such that \f$ W(r) = 0 \f$ if \f$ \require{AMSmath} \left\lvert r \right\rvert > h \f$, where \f$ W \f$ is the [kernel function](#kernel).

| Keyword                   | Description                                   | Unit             | Modes                        |
| ------------------------- | --------------------------------------------- | ---------------- | :--------------------------: |
| `smooth_names`            | Names of the smooth particles                 | string, no space | [&Dagger;](#Simulation mode) |
| `smooth_masses`           | Masses of the smooth particles                | kg               | [&Dagger;](#Simulation mode) |
| `smooth_unitmasses`       | Masses of a microscopic particle of the fluid | kg               | [&Dagger;](#Simulation mode) |
| `smooth_bulk_viscosity`   | Bulk viscosities \f$ \eta \f$ in the fluid    | Pa.s             | [&Dagger;](#Simulation mode) |
| `smooth_shear_viscosity`  | Shear viscosities \f$ \zeta \f$ in the fluid  | Pa.s             | [&Dagger;](#Simulation mode) |
| `smooth_smoothing_length` | Smoothing lengths \f$ h \f$                   | m                | [&Dagger;](#Simulation mode) |
| `smooth_kernel`           | Smoothing kernels \f$ W \f$                   | [Doc](#kernel)   | [&Dagger;](#Simulation mode) |

### <a name="kernel"></a> Kernel functions

Kernel functions are a key component of the [SDPD](#sdpd) method and allow to evaluate the field variables as an overage over neighboring particles. They are set through the keyword `smooth_kernel` (or `wall_kernel` for virtual wall particles). An associated parameter is the smoothing length \f$ h \f$ that determined the distance at which the kernel \f$ W \f$ vanishes. The keyword `smooth_smoothing_length` (respectively `wall_smoothing_length`) allows to fix it for each type of [smooth particles](#smoothparticle) or [wall particles](#wallparticle).

| Keyword: `smooth_kernel` | Name           | Link         | 
| ------------------------ | -------------- | :----------: |
| `lucy`                   | Lucy function. | [Doc](#lucy) |


#### <a name="lucy"></a> Lucy function

The Lucy function @cite lucy_1977 is a \f$ 4^{\rm th} \f$ order polynomial which reads
\f[
W_{\rm Lucy}(r) = \frac{105}{16\pi h^3}\left(1+3\frac{r}{h}\right)\left(1-\frac{r}{h}\right)^3 1_{0\leq r \leq h}.
\f]
For each smooth particle using it, the keyword `smooth_kernel` must be set to `lucy`.


## <a name="interaction"></a> 4. Interactions 

The SDPD particles usually do not interact through a potential function. It may however be useful to add one, for instance to model repulsive walls. We refer to [the MD help page](md_markdown_md.html#potential) for instructions on how to add a potential interaction.

### <a name="sdpd-interaction"></a> SDPD 

In order to activate both the conservative and the fluctuation/dissipation parts of SDPD, the following keywords are required.

| Keyword        | Description            | Unit             | 
| -------------- | ---------------------- | ---------------- | 
| `sdpd_type_A`  | Names of the particles | string, no space | 
| `sdpd_type_B`  | Names of the particles | string, no space | 

### <a name="wall-interaction"></a> Virtual SDPD for walls

This activates the computation of density and the SDPD conservative forces between real smooth particles (set with `smooth_types` keyword) and virtual wall particles (set with `wall_types`). The fluctuation/dissipation part of SDPD is not included. In order to ensure that real particles do not penetrate the wall, an additional repulsive potential may be added independently (for instance a [Lennard-Jones](#lennard-jones) potential). The handling of the virtual particles is described in @cite bian_2012 (and briefly in @cite faure_2016).

| Keyword                  | Description                    | Unit             | 
| ------------------------ | ------------------------------ | ---------------- | 
| `sdpd_wall_type_real`    | Names of the real particles    | string, no space | 
| `sdpd_wall_type_virtual` | Names of the virtual particles | string, no space | 


## <a name="eos"></a> 5. Equations of state

The thermodynamic state (entropy \f$ S_i \f$, temperature \f$ T_i \f$, pressure \f$ P_i \f$ and heat capacity \f$ C_i \f$) of a fluid particle is determined from the density and the internal energy thanks to an [equation of state](#eos). The following equation of states are compatible with SDPD. They are all given as the entropy function of the internal energy and the density. The parameters are specified thanks to array-valued input field summarized in the tables below.

### <a name="eos-ideal"></a> Ideal gas

The equation of state of an ideal gas is given by
\f[
\mathcal{S}(\varepsilon,\rho) = \frac{3}{2}Kk_{\rm B}{\rm ln}(\varepsilon) - \frac{1}{2}Kk_{\rm B}{\rm ln}(\rho).
\f]
The parameter \f$ K \f$ is the size of the smooth particle computed as the ratio between the mass \f$ m_i \f$ of the smooth particle (fixed with the keyword `smooth_masses`) and the mass of a microscopic particle of the fluid (fixed with the keyword `smooth_unitmasses`).

| Keyword   | Description                                                | Unit             |
| --------- | ---------------------------------------------------------- | ---------------- |
| `ig_type` | Name of the particles with the ideal gas equation of state | string, no space |

### <a name="mie-gruneisen"></a> Mie-Grüneisen

The Mie-Grüneisen equation of state @cite heuze_2012 is defined as
\f[
\mathcal{S}(\varepsilon,\rho) = C_{V_r}\log\left(1 + \frac{\varepsilon - \mathfrak{E}_k(\rho)}{C_{V_r}\Theta(\rho)u_r}\right),
\f]
with the reference energy
\f[
\mathfrak{E}_k(\rho) = \frac{K_s}{\rho_s}\left[\frac{\exp\left((N_s+1)\left[1-\frac{\rho_s}{\rho}\right]\right)}{(N_s+1)^2}-\frac{1-\frac{\rho_s}{\rho}}{N_s+1}\right],
\f]
and temperature
\f[
\Theta(\rho) = \Theta_0\left(\frac{\rho}{\rho_0}\right)^{\Gamma_{\infty}}\exp\left(\frac{\Gamma_0-\Gamma_{\infty}}{q}\left[1-\left(\frac{\rho_0}{\rho}\right)^q\right]\right).
\f]

| Keyword   | Description                                                    | Unit                         |
| --------- | -------------------------------------------------------------- | ---------------------------- |
| `mg_type` | Name of the particles with the Mie-Grüneisen equation of state | string, no space             |
| `mg_g0`   | Grüneisen parameter \f$ \gamma_0 \f$                           | -                            |
| `mg_ginf` | Grüneisen parameter \f$ \gamma_{\infty} \f$                    | -                            |
| `mg_t0`   | Reference temperature \f$ \Theta_0 \f$                         | K                            |
| `mg_q`    | Parameter \f$ q \f$                                            | -                            |
| `mg_rho0` | Reference density \f$ \rho_0 \f$                               | kg.m\f$^{-3}\f$              |
| `mg_ks`   | Parameter \f$ K_s \f$                                          | Pa                           |
| `mg_ns`   | Parameter \f$ n_s \f$                                          | -                            |
| `mg_rhos` | Reference density \f$ \rho_s \f$                               | kg.m\f$^{-3}\f$              |
| `mg_ur`   | Reference reduced energy \f$ u_r \f$                           | -                            |
| `mg_cvr`  | specific heat capacity \f$ C_{V_r} \f$                         | J.K\f$^{-1}\f$.kg\f$^{-1}\f$ |

### <a name="hz"></a> HZ

The HZ equationf of state is given by
\f[
\mathcal{S}(\varepsilon,\rho) = C_{V}\log\left[\frac{\varepsilon-\mathcal{E}_{\rm ref}(\rho)}{C_{V}} + \theta(\rho)\right] + C_{V}\Gamma_0\frac{\rho_0}{\rho},
\f]
with
\f[
\theta(\rho) = (T_0 - T_{00})\exp\left[\Gamma_0\left(1-\frac{\rho_0}{\rho}\right)\right],
\f]
and
\f[
\mathcal{E}_{\rm ref}(\rho) = \frac12\frac{c_0^2x^2}{1-sx} \times \left\{
\begin{array}{cl}
\displaystyle 1+\frac{sx}3-s\left(\Gamma_0-s\right)\frac{x^2}6 & \displaystyle \text{ if } x \geq 0,\\
\displaystyle 1 & \displaystyle \text{ if } x < 0,
\end{array}
\right.
\f]
where \f$ x=1-\frac{\rho_0}{\rho} \f$, and \f$T_0\f$ and \f$T_{00}\f$ are two constants defined as the standard temperature \f$T_0 = 298.13\f$ K and the temperature \f$T_{00}\f$ on the reference curve \f$\mathcal{E}_{\rm ref}\f$.
This constant is determined as \f$T_{00} = \frac{E_0}{C_V}\f$ where \f$E_0\f$ is the energy in standard conditions (temperature \f$T_0\f$ and pressure \f$P_0 = 10^5\f$~Pa).

| Keyword   | Description                                         | Unit                         |
| --------- | --------------------------------------------------- | ---------------------------- |
| `hz_type` | Name of the particles with the HZ equation of state | string, no space             |
| `hz_g0`   | Grüneisen parameter \f$ \gamma_0 \f$                | -                            |
| `hz_rho0` | Reference density \f$ \rho_0 \f$                    | kg.m\f$^{-3}\f$              |
| `hz_c0`   | Sound velocity \f$ c_0 \f$                          | m.s\f$^{-1}\f$               |
| `hz_cv`   | Specific heat capacity \f$ C_{V} \f$                | J.K\f$^{-1}\f$.kg\f$^{-1}\f$ |
| `hz_s`    | Parameter \f$ s \f$                                 | -                            |

### <a name="jwl"></a> JWL

The Jones-wilkins-Lee (JWL) equation of state reads
\f[
\mathcal{S}(\varepsilon,\rho) = C_{V}\log\left[\frac{\varepsilon-\mathcal{E}_k(\rho)}{C_{V}}\right] - C_{V}\Gamma_0\log(\rho), 
\f]
with
\f[
\mathcal{E}_k(\rho) = \frac{a}{\rho_0R_1}\exp\left[-R_1\frac{\rho_0}{\rho}\right] + \frac{b}{\rho_0R_2}\exp\left[-R_2\frac{\rho_0}{\rho}\right] + \frac{\mathfrak{K}}{\rho_0\Gamma_0}\left(\frac{\rho_0}{\rho}\right)^{-\Gamma_0} + C_{\rm ek},
\f]
In order to define the constants $\mathfrak{K}$ and $C_{\rm ek}$, we first introduce
\f{equation}{
\left\{
\begin{aligned}
\rho_{\rm CJ} &= \rho_0 \frac{\rho_0D_{\rm CJ}^2}{\rho_0D_{\rm CJ}^2-P_{\rm CJ}},\\
E_{\rm CJ} &= E_0 + \frac12P_{\rm CJ}\left(\frac1{\rho_0}-\frac1{\rho_{\rm CJ}}\right),\\
P_{\rm k1CJ} &= a\exp\left(-R_1\frac{\rho_0}{\rho_{\rm CJ}}\right) + b\exp\left(-R_2\frac{\rho_0}{\rho_{\rm CJ}}\right).
\end{aligned}
\right.
\f}
Then,
\f[
  \mathfrak{K} = \left(P_{\rm CJ} - P_{K1{\rm CJ}} - \frac1mC_v\Gamma_0T_{\rm CJ}\rho_{\rm CJ}\right)\left(\frac{\rho_0}{\rho_{\rm CJ}}\right)^{\Gamma_0+1},
\f]
and
\f[
  C_{\rm ek} = E_{\rm CJ} -\frac{a}{\rho_0R_1}\exp\left(-R_1\frac{\rho_0}{\rho_{\rm CJ}}\right) -\frac{b}{\rho_0R_2}\exp\left(-R_2\frac{\rho_0}{\rho_{\rm CJ}}\right) -\frac{P_{\rm CJ}-P_{\rm k1CJ}}{\rho_{\rm CJ}\Gamma_0}.
\f]

| Keyword    | Description                                          | Unit                         |
| ---------- | ---------------------------------------------------- | ---------------------------- |
| `jwl_type` | Name of the particles with the JWL equation of state | string, no space             |
| `jwl_g0`   | Grüneisen parameter \f$ \gamma_0 \f$                 | -                            |
| `jwl_rho0` | Reference density \f$ \rho_0 \f$                     | kg.m\f$^{-3}\f$              |
| `jwl_e0`   | Reference energy \f$E_0\f$                           | J                            |
| `jwl_dcj`  | Detonation velocity \f$D_{\rm CJ}\f$                 | m.s\f$^{-1}\f$               |
| `jwl_pcj`  | Pressure at CJ point \f$P_{\rm CJ}\f$                | Pa                           |
| `jwl_tcj`  | Temperature at CJ point \f$T_{\rm CJ}\f$             | K                            |
| `jwal_cv`  | specific heat capacity \f$C_{V}\f$                   | J.K\f$^{-1}\f$.kg\f$^{-1}\f$ |
| `jwl_a`    | Parameter \f$a\f$                                    | Pa                           |
| `jwl_b`    | Parameter \f$b\f$                                    | Pa                           |
| `jwl_r1`   | Parameter \f$r_1\f$                                  | -                            |
| `jwl_r2`   | Parameter \f$r_2\f$                                  | -                            |


## <a name="chemistry"></a> 6. Chemical reactions

In order to extend SDPD to reactive materials @cite faure_2018, an additional variable \f$ \lambda_i \in [0,1] \f$ is associated with each particle to describe the progress of the chemical reaction: \f$ A \to B \f$. When \f$ \lambda_i = 0 \f$, the mesoparticle is made entirely of the reactant \f$ A \f$ while \f$ \lambda_i=1 \f$ corresponds to a pure particle of the product \f$ B \f$. The progress variable \f$ \lambda_i \f$ represents the portion of the mesoparticle that has reacted and transformed from the reactant \f$ A \f$ into the product \f$ B\f$. 
The reactive model relies on two main ingredients: a [kinetics](#kinetics) which governs the evolution of the progress variable and a [reactive equation of state](#reactive) allowing to switch the equation of state from the product to the reactant along the simulation.

### <a name="kinetics"></a> Kinetics

The kinetics of the chemical reaction is responsible for the evolution of the progress variable. The exothermicity is also reinjected in the internal energy.

#### <a name="kinetics-so"></a> Second order kinetics

when a second order kinetics is considered, the progress variable is evolved through interaction with the neighboring particles as:
\f[
  \frac{d \lambda_i}{d t} = \sum_{j\neq i} \mathcal{K}_{0\to1}\left(T_{ij}\right)(1-\lambda_i)(1-\lambda_j)W(r_{ij}) - \mathcal{K}_{1\to0}\left(T_{ij}\right)\lambda_i\lambda_jW(r_{ij}),
\f]
where \f$ W \f$ is the [kernel function](#kernel) and \f$\mathcal{K}_{0\to1}\f$ and \f$\mathcal{K}_{1\to0}\f$ are the reaction rates, respectively, for the forward and backward reactions.
The reaction rates depend on the mean temperature \f$T_{ij} = \frac12\left(T_i+T_j\right)\f$ according to some Arrhenius law :
\f[
\mathcal{K}_{\mathcal{X}}(T) = Z_{\mathcal{X}}\exp\left(-\frac{E_{\mathcal{X}}}{k_{\rm B}T}\right),
\f]
with an activation energy \f$ E_{\mathcal{X}}\f$, that represents the energy barrier a molecule needs to overcome during the reaction, and a prefactor \f$Z_{\mathcal{X}}\f$ that governs the frequency of the reaction.

The kinetics is discretized with an explicit Euler scheme and the exothermicity \f$ E_{\rm exo} = E_{1\to0}- E_{0\to1} \f$ is reinjected in the internal energy so as to preserve the total energy now reading
\f[
  E(q,p,\varepsilon,\lambda) = \sum_{i=1}^N \varepsilon_i + \frac{p_i^2}{2m} + (1-\lambda_i)KE_{\rm exo}.
\f]
Using a time step \f$ \Delta t \f$, the updated progress variables and internal energies are given by
\f{equation}{
\left\{
\begin{aligned}
\lambda_i^{n+1} &= \lambda_i^n + \sum_{j\neq i} \mathcal{K}_{0\to1}\left(T_{ij}^n\right)(1-\lambda_i^n)(1-\lambda_j^n)W(r_{ij}^n) - \mathcal{K}_{1\to0}\left(T_{ij}^n\right)\lambda_i^n\lambda_j^nW(r_{ij}^n),\\
\varepsilon_i^{n+1} &= \varepsilon_i^n + (\lambda_i^{n+1}-\lambda_i^n)KE_{\rm exo}.
\end{aligned}
\right.
\f}

| Keyword                           | Description                                                   | Unit             | 
| --------------------------------- | ------------------------------------------------------------- | ---------------- | 
| `reaction_so_type`                | Name of the particles subject to a chemical reaction          | string, no space | 
| `reaction_so_forward_prefactor`   | Prefactor for the forward reaction \f$ Z_{0\to1} \f$          | s$f^{-1}\f$      | 
| `reaction_so_backward_prefactor`  | Prefactor for the backward reaction \f$ Z_{0\to1} \f$         | s$f^{-1}\f$      | 
| `reaction_so_forward_activation`  | Activation energy for the forward reaction \f$ Z_{0\to1} \f$  | J                |
| `reaction_so_backward_activation` | Activation energy for the backward reaction \f$ Z_{0\to1} \f$ | J                |


### <a name="reactive"></a> Reactive equation of state

The reactive equation of state allows to switch the equation of state for a type of particle along the simulation. An equation of state for both the reactant and the product must be defined independently among the usual [equations of state](#eos). The internal energy and the density of a mixed mesoparticle are given by
\f{equation}{
\begin{aligned}
\varepsilon_i &= \varepsilon_i^0 + \varepsilon_j^0,\\
\rho_i &= (1-\lambda)\rho_i^0 + \lambda\rho_i^1,
\end{aligned}
\f}
where the superscript denotes the quantities associated with the reactant (\f$ 0 \f$) or the product (\f$ 1 \f$).
The reactive equation of state mixes the two equations of state so that thermal and mechanical equilibrium is achieved inside a mesoparticle, *i.e*
\f{equation}{
\begin{aligned}
\mathcal{T}^0(\varepsilon_i^0,\rho_i^0) &= \mathcal{T}^1(\varepsilon_i^1,\rho_i^1), \\
\mathcal{P}^0(\varepsilon_i^0,\rho_i^0) &= \mathcal{P}^1(\varepsilon_i^1,\rho_i^1).
\end{aligned}
\f}
A Newton algorithm is used to perform the numerical inversion.

| Keyword         | Description                                          | Unit             | 
| --------------- | ---------------------------------------------------- | ---------------- | 
| `reactive_type` | Name of the particles subject to a chemical reaction | string, no space | 
| `reactive_eos0` | Name of the equation of state of the reactant        | string           | 
| `reactive_eos1` | Name of the equation of state of the product         | string           | 

The names of the equation of state for the reactant and product are to be chosen in the following table:

| Keyword: `reactive_eos0`/`reactive_eos1` | Description                                                                    | 
| ---------------------------------------- | ------------------------------------------------------------------------------ |
| `HZ`                                     | [HZ](#hz) equation of state                                                    |
| `Jones-Wilkins-Lee`                      | [JWL](#jwl) equation of state                                                  |

## <a name="domain"></a> 7. Setting the domain

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

#### <a name="wall"></a> Walls

Walls are modelled as a repulsive Lennard-Jones potential with \f$ \sigma \approx R_{\rm cut}\f$, with \f$ R_{\rm cut} \f$ the maximum cut-off radius in the simulation (potential or kernel smoothing length). However it is also possible to model the walls as a lattice of fixed particles. The interaction with the other (real) particles then needs to be specified as any other particle types. It may be with a [potential](md_markdown_md.html#potential) or with [virtual SDPD interactions](#wall-interaction).

When initializing with a [lattice](#lattice), the width of the layer is specified thanks to the keyword `wall_width`. All lattice sites which are closer to the domain limit will be populated with wall particles (if the `wall` boundary conditions has been chosen for this direction). For each direction, the type of wall particles to be used is chosen with `wall_lower_types` for the lower end of the domain and `wall_upper_types` for the upper end). The value `null` means that no particle is put within the wall.

| Keyword                 | Description                                                                   | Unit                 |
| ----------------------- | ----------------------------------------------------------------------------- | -------------------- |
| `wall_names`            | Names of the wall particles                                                   | string, no space     |
| `wall_masses`           | Masses of the wall particles                                                  | kg                   |
| `wall_unitmasses`       | Masses of a microscopic particle in the wall                                  | kg                   |
| `wall_smoothing_length` | Smoothing lengths                                                             | m                    |
| `wall_kernel`           | Kernel functions                                                              | string, no space     |
| `wall_velocity_x`       | Constant velocities in the x-direction                                        | m.s\f$^{-1}\f$       |
| `wall_velocity_y`       | Constant velocities in the x-direction                                        | m.s\f$^{-1}\f$       |
| `wall_velocity_z`       | Constant velocities in the x-direction                                        | m.s\f$^{-1}\f$       |
| `wall_width`            | Width of the wall                                                             | m                    |
| `wall_lower_types`      | Types of wall particles for each direction at the lower end of the domain     | string (3-dim array) |
| `wall_upper_types`      | Types of wall particles for each direction at the upper end of the domain     | string (3-dim array) |
| `wall_stop_time`        | Time limit after which the walls are stopped (-1 to let them go indefinitely) | s                    |

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

## <a name="output"></a> 8. Output

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

## <a name="parallel"></a> 9. Parallelisation

| Keyword                | Description                               |
| ---------------------- | ----------------------------------------- |
| `decoupage`            | Number of MPI processes in each direction |
| `max_threads_per_node` | Maximum number of TBB threads per process |

# <a name="example"></a> Examples

## Input file for SDPD
@include ../../tests/single-mat-sdpd.xsp

## Input file for SDPD with walls (shock)
@include ../../tests/single-mat-sdpd-shock.xsp

## Input file for reactive SDPD
@include ../../tests/single-mat-sdpd-reactive.xsp



