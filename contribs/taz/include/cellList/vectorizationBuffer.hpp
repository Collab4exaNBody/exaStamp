/// @file
/// @brief Definition of the class VectBuffer

#ifndef __VECTORIZATION_BUFFER_HPP_INCLUDED
#define __VECTORIZATION_BUFFER_HPP_INCLUDED


#include "utils/array/extArray.hpp"


/// @brief Provide memory-aligned arrays to perform vectorized operations in
/// force computation and neighbor list build
/// @tparam align Alignment
template <size_t align> class VectBuffer {
public:

  /// @brief Default constructor
  VectBuffer() {}

  /// @brief Destructor (nothing to do)
  ~VectBuffer() {}

  /// @brief Get a pointer to the particle indices
  inline uint64_t* id() { return m_id.data(); }

  /// @brief Get a pointer to the x components of the distances
  inline double* drx() { return m_drx.data(); }

  /// @brief Get a pointer to the y components of the distances
  inline double* dry() { return m_dry.data(); }

  /// @brief Get a pointer to the z components of the distances
  inline double* drz() { return m_drz.data(); }

  /// @brief Get a pointer to the x components of the velocities
  inline double* dvx() { return m_dvx.data(); }
  /// @brief Get a pointer to the y components of the velocities
  inline double* dvy() { return m_dvy.data(); }
  /// @brief Get a pointer to the z components of the velocities
  inline double* dvz() { return m_dvz.data(); }

  /// @brief Get a pointer to the screening term
  inline double* S() { return m_S.data(); }

  /// @brief Get a pointer to new index of neighbaros
  inline uint* indxBuffer() { return m_indxBuffer.data(); }

  /// @brief Get a pointer to the derivative of screening term
  inline double* dfSx() { return m_dfSx.data(); }

  /// @brief Get a pointer to the derivative of screening term
  inline double* dfSy() { return m_dfSy.data(); }

  /// @brief Get a pointer to the derivative of screening term
  inline double* dfSz() { return m_dfSz.data(); }

  /// @brief Get a pointer to the derivative of rho0
  inline double* drho0x() { return m_drho0x.data(); }

  /// @brief Get a pointer to the derivative of rho0
  inline double* drho0y() { return m_drho0y.data(); }

  /// @brief Get a pointer to the derivative of rho0
  inline double* drho0z() { return m_drho0z.data(); }

  /// @brief Get a pointer to the derivative of rho1
  inline double* drho1x() { return m_drho1x.data(); }

  /// @brief Get a pointer to the derivative of rho1
  inline double* drho1y() { return m_drho1y.data(); }

  /// @brief Get a pointer to the derivative of rho1
  inline double* drho1z() { return m_drho1z.data(); }

  /// @brief Get a pointer to the derivative of rho2
  inline double* drho2x() { return m_drho2x.data(); }

  /// @brief Get a pointer to the derivative of rho2
  inline double* drho2y() { return m_drho2y.data(); }

  /// @brief Get a pointer to the derivative of rho2
  inline double* drho2z() { return m_drho2z.data(); }

  /// @brief Get a pointer to the derivative of rho3
  inline double* drho3x() { return m_drho3x.data(); }

  /// @brief Get a pointer to the derivative of rho3
  inline double* drho3y() { return m_drho3y.data(); }

  /// @brief Get a pointer to the derivative of rho3
  inline double* drho3z() { return m_drho3z.data(); }

  /// @brief Get a pointer to the x components of the forces
  inline double* dfx() { return m_dfx.data(); }

  /// @brief Get a pointer to the y components of the forces
  inline double* dfy() { return m_dfy.data(); }

  /// @brief Get a pointer to the z components of the forces
  inline double* dfz() { return m_dfz.data(); }

  /// @brief Get a pointer to the x potential energies
  inline double* den() { return m_den.data(); }

  /// @brief Get a pointer to the rho terms
  inline double* rho() { return m_tmp0.data(); }

  /// @brief Get a pointer to the embedding terms
  inline double* emb() { return m_tmp1.data(); }

  /// @brief Get a pointer to the rcut data
  inline double* rct() { return m_tmp0.data(); }

  /// @brief Get a pointer to the neighbor computation results
  inline double* out() { return m_tmp1.data(); }

  /// @brief Get a pointer to the local inverse temperature \f$\beta_{i,j}=k_B\left(\frac1{T_i}+\frac1{T_j}\right) \f$
  inline double* bij() { return m_tmp2.data(); }
  /// @brief Get a pointer to the SDPD fluctuation amplitude
  inline double* sij() { return m_tmp3.data(); }


  /// @brief Get a constant pointer to the particle indices
  inline const uint64_t* id() const { return m_id.data(); }

  /// @brief Get a constant pointer to the x components of the distances
  inline const double* drx() const { return m_drx.data(); }

  /// @brief Get a constant pointer to the y components of the distances
  inline const double* dry() const { return m_dry.data(); }

  /// @brief Get a constant pointer to the z components of the distances
  inline const double* drz() const { return m_drz.data(); }

  /// @brief Get a constant pointer to the x components of the velocities
  inline const double* dvx() const { return m_dvx.data(); }
  /// @brief Get a constant pointer to the x components of the velocities
  inline const double* dvy() const { return m_dvy.data(); }
  /// @brief Get a constant pointer to the x components of the velocities
  inline const double* dvz() const { return m_dvz.data(); }

  /// @brief Get a constant pointer to the screening term
  inline const double* S() const { return m_S.data(); }

  /// @brief Get a pointer to new index of neighbaros
  inline const uint* indxBuffer() const { return m_indxBuffer.data(); }

  /// @brief Get a constant pointer to the screening term
  inline const double* dfSx() const { return m_dfSx.data(); }

  /// @brief Get a constant pointer to the screening term
  inline const double* dfSy() const { return m_dfSy.data(); }

  /// @brief Get a constant pointer to the screening term
  inline const double* dfSz() const { return m_dfSz.data(); }

  /// @brief Get a constant pointer to rho0
  inline const double* drho0x() const { return m_drho0x.data(); }

  /// @brief Get a constant pointer to rho0
  inline const double* drho0y() const { return m_drho0y.data(); }

  /// @brief Get a constant pointer to rho0
  inline const double* drho0z() const { return m_drho0z.data(); }

  /// @brief Get a constant pointer to rho1
  inline const double* drho1x() const { return m_drho1x.data(); }

  /// @brief Get a constant pointer to rho1
  inline const double* drho1y() const { return m_drho1y.data(); }

  /// @brief Get a constant pointer to rho1
  inline const double* drho1z() const { return m_drho1z.data(); }

  /// @brief Get a constant pointer to rho2
  inline const double* drho2x() const { return m_drho2x.data(); }

  /// @brief Get a constant pointer to rho2
  inline const double* drho2y() const { return m_drho2y.data(); }

  /// @brief Get a constant pointer to rho2
  inline const double* drho2z() const { return m_drho2z.data(); }

  /// @brief Get a constant pointer to rho3
  inline const double* drho3x() const { return m_drho3x.data(); }

  /// @brief Get a constant pointer to rho3
  inline const double* drho3y() const { return m_drho3y.data(); }

  /// @brief Get a constant pointer to rho3
  inline const double* drho3z() const { return m_drho3z.data(); }

  /// @brief Get a constant pointer to the x components of the forces
  inline const double* dfx() const { return m_dfx.data(); }

  /// @brief Get a constant pointer to the y components of the forces
  inline const double* dfy() const { return m_dfy.data(); }

  /// @brief Get a constant pointer to the z components of the forces
  inline const double* dfz() const { return m_dfz.data(); }

  /// @brief Get a constant pointer to the potential energies
  inline const double* den() const { return m_den.data(); }

  /// @brief Get a constant pointer to the rho terms
  inline const double* rho() const { return m_tmp0.data(); }

  /// @brief Get a constant pointer to the embedding terms
  inline const double* emb() const { return m_tmp1.data(); }

  /// @brief Get a constant pointer to the rcut data
  inline const double* rct() const { return m_tmp0.data(); }

  /// @brief Get a constant pointer to the neighbor computation results
  inline const double* out() const { return m_tmp1.data(); }

  /// @brief Get a constant pointer to the local inverse temperature \f$\beta_{i,j}=k_B\left(\frac1{T_i}+\frac1{T_j}\right) \f$
  inline const double* bij() const { return m_tmp2.data(); }
  /// @brief Get a constant pointer to the SDPD fluctuation amplitude
  inline const double* sij() const { return m_tmp3.data(); }

  
  /// @brief Get a reference to the global index of the distance for a specified particle
  /// @param [in] i Particle index
  inline uint64_t& id(const uint i) { return m_id[i]; }

  /// @brief Get a reference to screening term
  /// @param [in] i Particle index
  inline double& S(const uint i) { return m_S[i]; }

  /// @brief Get a pointer to new index of neighbours
  inline uint& indxBuffer(const uint i) { return m_indxBuffer[i]; }

  /// @brief Get a reference to the derivative of screening term
  /// @param [in] i Particle index
  inline double& dfSx(const uint i) { return m_dfSx[i]; }

  /// @brief Get a reference to the derivative of screening term
  /// @param [in] i Particle index
  inline double& dfSy(const uint i) { return m_dfSy[i]; }

  /// @brief Get a reference to the derivative of screening term
  /// @param [in] i Particle index
  inline double& dfSz(const uint i) { return m_dfSz[i]; }

  /// @brief Get a reference to the derivative of rho0
  /// @param [in] i Particle index
  inline double& drho0x(const uint i) { return m_drho0x[i]; }

  /// @brief Get a reference to the derivative of rho0
  /// @param [in] i Particle index
  inline double& drho0y(const uint i) { return m_drho0y[i]; }

  /// @brief Get a reference to the derivative of rho0
  /// @param [in] i Particle index
  inline double& drho0z(const uint i) { return m_drho0z[i]; }

  /// @brief Get a reference to the derivative of rho1
  /// @param [in] i Particle index
  inline double& drho1x(const uint i) { return m_drho1x[i]; }

  /// @brief Get a reference to the derivative of rho1
  /// @param [in] i Particle index
  inline double& drho1y(const uint i) { return m_drho1y[i]; }

  /// @brief Get a reference to the derivative of rho1
  /// @param [in] i Particle index
  inline double& drho1z(const uint i) { return m_drho1z[i]; }

  /// @brief Get a reference to the derivative of rho2
  /// @param [in] i Particle index
  inline double& drho2x(const uint i) { return m_drho2x[i]; }

  /// @brief Get a reference to the derivative of rho2
  /// @param [in] i Particle index
  inline double& drho2y(const uint i) { return m_drho2y[i]; }

  /// @brief Get a reference to the derivative of rho2
  /// @param [in] i Particle index
  inline double& drho2z(const uint i) { return m_drho2z[i]; }

  /// @brief Get a reference to the derivative of rho3
  /// @param [in] i Particle index
  inline double& drho3x(const uint i) { return m_drho3x[i]; }

  /// @brief Get a reference to the derivative of rho3
  /// @param [in] i Particle index
  inline double& drho3y(const uint i) { return m_drho3y[i]; }

  /// @brief Get a reference to the derivative of rho3
  /// @param [in] i Particle index
  inline double& drho3z(const uint i) { return m_drho3z[i]; }

  /// @brief Get a reference to the x component of the distance for a specified particle
  /// @param [in] i Particle index
  inline double& drx(const uint i) { return m_drx[i]; }

  /// @brief Get a reference to the y component of the distance for a specified particle
  /// @param [in] i Particle index
  inline double& dry(const uint i) { return m_dry[i]; }

  /// @brief Get a reference to the z component of the distance for a specified particle
  /// @param [in] i Particle index
  inline double& drz(const uint i) { return m_drz[i]; }

  /// @brief Get a reference to the x component of the velocity for a specified particle
  /// @param [in] i Particle index
  inline double& dvx(const uint i) { return m_dvx[i]; }
  /// @brief Get a reference to the y component of the velocity for a specified particle
  /// @param [in] i Particle index
  inline double& dvy(const uint i) { return m_dvy[i]; }
  /// @brief Get a reference to the z component of the velocity for a specified particle
  /// @param [in] i Particle index
  inline double& dvz(const uint i) { return m_dvz[i]; }

  /// @brief Get a reference to the x component of the force for a specified particle
  /// @param [in] i Particle index
  inline double& dfx(const uint i) { return m_dfx[i]; }

  /// @brief Get a reference to the y component of the force for a specified particle
  /// @param [in] i Particle index
  inline double& dfy(const uint i) { return m_dfy[i]; }

  /// @brief Get a reference to the z component of the force for a specified particle
  /// @param [in] i Particle index
  inline double& dfz(const uint i) { return m_dfz[i]; }

  /// @brief Get a reference to the potential energy for a specified particle
  /// @param [in] i Particle index
  inline double& den(const uint i) { return m_den[i]; }

  /// @brief Get a reference to the rho term for a specified particle
  /// @param [in] i Particle index
  inline double& rho(const uint i) { return m_tmp0[i]; }

  /// @brief Get a reference to the embedding term for a specified particle
  /// @param [in] i Particle index
  inline double& emb(const uint i) { return m_tmp1[i]; }

  /// @brief Get a reference to the rcut for a specified particle
  /// @param [in] i Particle index
  inline double& rct(const uint i) { return m_tmp0[i]; }

  /// @brief Get a reference to the neighbor computation result for a specified particle
  /// @param [in] i Particle index
  inline double& out(const uint i) { return m_tmp1[i]; }

  /// @brief Get a reference to the local inverse temperature for a specified particle
  inline double& bij(const uint i) { return m_tmp2[i]; }
  /// @brief Get a reference to the SDPD fluctuation amplitude for a specified particle
  inline double& sij(const uint i) { return m_tmp3[i]; }



  /// @brief Get a constant reference to the global index for a specified particle
  /// @param [in] i Particle index
  inline const uint64_t& id(const uint i) const { return m_id[i]; }

  /// @brief Get a constant reference to the x component of the distance for a specified particle
  /// @param [in] i Particle index
  inline const double& drx(const uint i) const { return m_drx[i]; }

  /// @brief Get a constant reference to the y component of the distance for a specified particle
  /// @param [in] i Particle index
  inline const double& dry(const uint i) const { return m_dry[i]; }

  /// @brief Get a constant reference to the z component of the distance for a specified particle
  /// @param [in] i Particle index
  inline const double& drz(const uint i) const { return m_drz[i]; }

  /// @brief Get a constant reference to the x component of the velocity for a specified particle
  /// @param [in] i Particle index
  inline const double& dvx(const uint i) const { return m_dvx[i]; }
  /// @brief Get a constant reference to the y component of the velocity for a specified particle
  /// @param [in] i Particle index
  inline const double& dvy(const uint i) const { return m_dvy[i]; }
  /// @brief Get a constant reference to the z component of the velocity for a specified particle
  /// @param [in] i Particle index
  inline const double& dvz(const uint i) const { return m_dvz[i]; }

  /// @brief Get a constant reference to the x component of the force for a specified particle
  /// @param [in] i Particle index
  inline const double& dfx(const uint i) const { return m_dfx[i]; }

  /// @brief Get a constant reference to the y component of the force for a specified particle
  /// @param [in] i Particle index
  inline const double& dfy(const uint i) const { return m_dfy[i]; }

  /// @brief Get a constant reference to the z component of the force for a specified particle
  /// @param [in] i Particle index
  inline const double& dfz(const uint i) const { return m_dfz[i]; }

  /// @brief Get a constant reference to the potential energy for a specified particle
  /// @param [in] i Particle index
  inline const double& den(const uint i) const { return m_den[i]; }

  /// @brief Get a constant reference to the rho term for a specified particle
  /// @param [in] i Particle index
  inline const double& rho(const uint i) const { return m_tmp0[i]; }

  /// @brief Get a constant reference to the embedding term for a specified particle
  /// @param [in] i Particle index
  inline const double& emb(const uint i) const { return m_tmp1[i]; }

  /// @brief Get a constant reference to the rcut for a specified particle
  /// @param [in] i Particle index
  inline const double& rct(const uint i) const { return m_tmp0[i]; }

  /// @brief Get a constant reference to the neighbor computation result for a specified particle
  /// @param [in] i Particle index
  inline const double& out(const uint i) const { return m_tmp1[i]; }

  /// @brief Get a constant reference to the local inverse temperature for a specified particle
  inline const double& bij(const uint i) const { return m_tmp2[i]; }
  /// @brief Get a constant reference to the SDPD fluctuation amplitude for a specified particle
  inline const double& sij(const uint i) const { return m_tmp3[i]; }

  /// @brief Used to stock screening term between particles j and k in the stockage of particle i.
  /// @param [in] i Particle index
  inline const double& S(const uint i) const { return m_S[i]; }

  /// @brief Get a pointer to new index of neighbaros
  inline const uint& indxBuffer(const uint i) const { return m_indxBuffer[i]; }

  /// @brief Used to stock partial contribution of screening term between particles j and k in the stockage of particle i.
  /// @param [in] i Particle index
  inline const double& dfSx(const uint i) const { return m_dfSx[i]; }

  /// @brief Used to stock partial contribution of screening term between particles j and k in the stockage of particle i.
  /// @param [in] i Particle index
  inline const double& dfSy(const uint i) const { return m_dfSy[i]; }

  /// @brief Used to stock partial contribution of screening term between particles j and k in the stockage of particle i.
  /// @param [in] i Particle index
  inline const double& dfSz(const uint i) const { return m_dfSz[i]; }

  /// @brief Used to stock rho1 
  /// @param [in] i Particle index
  inline const double& drho0x(const uint i) const { return m_drho0x[i]; }

  /// @brief Used to stock partial contribution of screening term between particles j and k in the stockage of particle i.
  /// @param [in] i Particle index
  inline const double& drho0y(const uint i) const { return m_drho0y[i]; }

  /// @brief Used to stock rho0
  /// @param [in] i Particle index
  inline const double& drho0z(const uint i) const { return m_drho0z[i]; }

  /// @brief Used to stock rho1 
  /// @param [in] i Particle index
  inline const double& drho1x(const uint i) const { return m_drho1x[i]; }

  /// @brief Used to stock partial contribution of screening term between particles j and k in the stockage of particle i.
  /// @param [in] i Particle index
  inline const double& drho1y(const uint i) const { return m_drho1y[i]; }

  /// @brief Used to stock rho1
  /// @param [in] i Particle index
  inline const double& drho1z(const uint i) const { return m_drho1z[i]; }

  /// @brief Used to stock rho2 
  /// @param [in] i Particle index
  inline const double& drho2x(const uint i) const { return m_drho2x[i]; }

  /// @brief Used to stock partial contribution of screening term between particles j and k in the stockage of particle i.
  /// @param [in] i Particle index
  inline const double& drho2y(const uint i) const { return m_drho2y[i]; }

  /// @brief Used to stock rho2
  /// @param [in] i Particle index
  inline const double& drho2z(const uint i) const { return m_drho2z[i]; }

  /// @brief Used to stock rho3 
  /// @param [in] i Particle index
  inline const double& drho3x(const uint i) const { return m_drho3x[i]; }

  /// @brief Used to stock partial contribution of screening term between particles j and k in the stockage of particle i.
  /// @param [in] i Particle index
  inline const double& drho3y(const uint i) const { return m_drho3y[i]; }

  /// @brief Used to stock rho3
  /// @param [in] i Particle index
  inline const double& drho3z(const uint i) const { return m_drho3z[i]; }



  /// @brief Resize all common arrays
  /// @param [in] n New size
  inline void resize(const uint n) {

    // Only if size increase
    if (n > m_drx.size()) {

      m_indxBuffer.resize(n); // for Verlet Lists
      m_id.resize(n);

      m_drx.resize(n);
      m_dry.resize(n);
      m_drz.resize(n);

      m_dvx.resize(n);
      m_dvy.resize(n);
      m_dvz.resize(n);
      
      m_dfx.resize(n);
      m_dfy.resize(n);
      m_dfz.resize(n);
      m_den.resize(n);

      m_tmp0.resize(n);
      m_tmp1.resize(n);
      m_tmp2.resize(n);
      m_tmp3.resize(n);

    }

  }


  /// @brief Resize all arrays used by MEAM Potential
  /// @param [in] n New size
  inline void resizeMeam(const uint n) {

    // Only if size increase
    if (n > m_dfSx.size()) {

      m_indxBuffer.resize(n);
      m_S.resize(n);
      m_dfSx.resize(n);
      m_dfSy.resize(n);
      m_dfSz.resize(n);

      m_drho0x.resize(n);
      m_drho0y.resize(n);
      m_drho0z.resize(n);

      m_drho1x.resize(n);
      m_drho1y.resize(n);
      m_drho1z.resize(n);

      m_drho2x.resize(n);
      m_drho2y.resize(n);
      m_drho2z.resize(n);

      m_drho3x.resize(n);
      m_drho3y.resize(n);
      m_drho3z.resize(n);

    	resize(n);

    }

  }

private:

  static constexpr uint8_t chunk = 16; ///< Chunk for the arrays

  std::vector<uint64_t> m_id; ///< Storage for the particles id

  std::vector<double> m_drx; ///< Storage for the x components of the distances
  std::vector<double> m_dry; ///< Storage for the y components of the distances
  std::vector<double> m_drz; ///< Storage for the z components of the distances

  std::vector<double> m_dvx; ///< Storage for the x components of the velocities
  std::vector<double> m_dvy; ///< Storage for the y components of the velocities
  std::vector<double> m_dvz; ///< Storage for the z components of the velocities

  std::vector<double> m_S; ///< Storage for the screening term
  std::vector<uint> m_indxBuffer; ///< Storage for new index of neighbors

  std::vector<double> m_dfSx; ///< Storage for the x components of derivative of screening term for neighborS
  std::vector<double> m_dfSy; ///< Storage for the y components of derivative of screening term for neighborS
  std::vector<double> m_dfSz; ///< Storage for the z components of derivative of screening term for neighborS

  std::vector<double> m_drho0x; ///< Storage for the x components of derivative of rho0
  std::vector<double> m_drho0y; ///< Storage for the y components of derivative of rho0
  std::vector<double> m_drho0z; ///< Storage for the z components of derivative of rho0

  std::vector<double> m_drho1x; ///< Storage for the x components of derivative of rho1
  std::vector<double> m_drho1y; ///< Storage for the y components of derivative of rho1
  std::vector<double> m_drho1z; ///< Storage for the z components of derivative of rho1

  std::vector<double> m_drho2x; ///< Storage for the x components of derivative of rho2
  std::vector<double> m_drho2y; ///< Storage for the y components of derivative of rho2
  std::vector<double> m_drho2z; ///< Storage for the z components of derivative of rho2

  std::vector<double> m_drho3x; ///< Storage for the x components of derivative of rho3
  std::vector<double> m_drho3y; ///< Storage for the y components of derivative of rho3
  std::vector<double> m_drho3z; ///< Storage for the z components of derivative of rho3

  std::vector<double> m_dfx; ///< Storage for the x components of the forces
  std::vector<double> m_dfy; ///< Storage for the y components of the forces
  std::vector<double> m_dfz; ///< Storage for the z components of the forces
  std::vector<double> m_den; ///< Storage for the potential energies

  std::vector<double> m_tmp0; ///< Storage for the rho term or rcut (depending on the calculation)
  std::vector<double> m_tmp1; ///< Storage for the embedding term or neighbor computation results (depending on the calculation)
  std::vector<double> m_tmp2; ///< Storage for the local inverse temperature
  std::vector<double> m_tmp3; ///< Storage for the SDPD fluctuation amplitudes

};


#endif // __VECTORIZATION_BUFFER_HPP_INCLUDED
