#pragma once

#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/particle_id_constants.h>
#include <onika/soatl/field_tuple.h>

namespace exaStamp
{
  using exaStamp::SlipLinkField;

    // rx,ry,rz : sliplink positon ( 3 reals )
    // vx,vy,vz : relative position of anchor to sliplink position
    // ep : fractional part of curvilinear position along the chain ( real )
    // cmol : [ left bead id, right bead id, sl id, peer sl id ] ( 4 integers )
    using SLGrid       =              Grid< FieldSet< /* rx,ry,rz implicitly included */  field::_id, field::_sl > >;
    using SLTupleInput = onika::soatl::FieldTuple< field::_rx, field::_ry, field::_rz, field::_id, field::_sl >;

    // constants defining the meaning of for integers stored in cmol
    static constexpr size_t SL_CMOL_LEFT_BEAD_ID = 0;
    static constexpr size_t SL_CMOL_RIGHT_BEAD_ID = 1;
    static constexpr size_t SL_CMOL_SL_ID = 2;
    static constexpr size_t SL_CMOL_SL_PEER_ID = 3;

    struct SlipLinkParameters
    {
      size_t number_of_chains = 1000;
      size_t beads_per_chain = 100;
      size_t number_of_sliplinks = 1000;
      double beta = 1.0;
      double Ns = 0.5;
      double kBT = 1.0;
      double ksi = 1.0;
      double sigma_rcut_scale = 5.0;
      
      double h;       // a.k.a dt (timestep)
      double ksis;

      double cte2;
      double cte4;
      double cte7;
      double cte8;
      double sigma1;
      double sigma2;
      double sigma3;
      double sigma4;
      double bond_max_dist;
    };

    struct SLBeadPlacement
    {
      uint64_t sl_id; // sliplink id
      uint64_t bead_id; // left bead id
      double frac; // fractional part of curvilinear coordinate, in [0;1[
    };

    struct SLBeadFrac
    {
      uint64_t bead_id; // left bead id
      double frac; // fractional part of curvilinear coordinate, in [0;1[     
    };

    struct alignas(8) SLRegenInfoNoId
    {
      uint64_t left_bead_id; // left bead id
      double xj_frac; // fractional part of curvilinear coordinate, in [0;1[
      bool full_regen;
    };

    struct alignas(8) SLRegenInfo
    {
      uint64_t sl_id;
      uint64_t left_bead_id; // left bead id
      double xj_frac; // fractional part of curvilinear coordinate, in [0;1[
      bool full_regen;
    };

    struct alignas(8) SLRegenInfo2
    {
      uint64_t bead_id; // left bead id
      Vec3d pos;
      double frac; // fractional part of curvilinear coordinate, in [0;1[
      bool full_regen;
    };

    static inline uint64_t sliplink_left_bead_id( uint64_t bead_id , size_t n_beads )
    {
      assert( bead_id != PARTICLE_MISSING_ID ); // only meaningful for local ids      
      if( bead_id==PARTICLE_NO_ID || ( (bead_id%n_beads) == 0 ) )
      {
        return PARTICLE_NO_ID;
      }
      else
      {
        return bead_id-1;
      }
    }

    static inline uint64_t sliplink_right_bead_id( uint64_t bead_id , size_t n_beads )
    {
      assert( bead_id != PARTICLE_MISSING_ID ); // only meaningful for local ids
      if( bead_id==PARTICLE_NO_ID || ( (bead_id%n_beads) == (n_beads-1) ) )
      {
        return PARTICLE_NO_ID;
      }
      else
      {
        return bead_id+1;
      }
    }

}

