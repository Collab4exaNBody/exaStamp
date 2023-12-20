#pragma once

#include <exanb/core/config.h>
#include <onika/memory/allocator.h>

namespace exaStamp
{

  using onika::memory::DEFAULT_ALIGNMENT;

  // additional storage space added to compute buffer created by compute_cell_particle_pairs
  struct alignas(DEFAULT_ALIGNMENT) RangeNeighborsExtraStorage
  {
    alignas(DEFAULT_ALIGNMENT) uint16_t nbh_list0[exanb::MAX_PARTICLE_NEIGHBORS];
    alignas(DEFAULT_ALIGNMENT) uint16_t nbh_list1[exanb::MAX_PARTICLE_NEIGHBORS];
    alignas(DEFAULT_ALIGNMENT) uint16_t nbh_list2[exanb::MAX_PARTICLE_NEIGHBORS];
    alignas(DEFAULT_ALIGNMENT) uint8_t nbh_type[exanb::MAX_PARTICLE_NEIGHBORS];
    uint16_t nbh_count0 = 0;
    uint16_t nbh_count1 = 0;
    uint16_t nbh_count2 = 0;

    inline unsigned int number_of_nbh0() const { return nbh_count0; }
    inline unsigned int number_of_nbh1() const { return nbh_count1; }
    inline unsigned int number_of_nbh2() const { return nbh_count2; }

    inline unsigned int nbh0(unsigned int i) const { return nbh_list0[i]; }
    inline unsigned int nbh1(unsigned int i) const { return nbh_list1[i]; }
    inline unsigned int nbh2(unsigned int i) const { return nbh_list2[i]; }
    
    inline bool nbh_is_type0(unsigned int i) const
    {
      const unsigned int shift = (i%2!=0) ? 4 : 0;
      return ( ( nbh_type[i/2] >> shift ) & 0x01 ) != 0 ;
    }

    inline bool nbh_is_type1(unsigned int i) const
    {
      const unsigned int shift = (i%2!=0) ? 4 : 0;
      return ( ( nbh_type[i/2] >> shift ) & 0x02 ) != 0 ;
    }

    inline bool nbh_is_type2(unsigned int i) const
    {
      const unsigned int shift = (i%2!=0) ? 4 : 0;
      return ( ( nbh_type[i/2] >> shift ) & 0x04 ) != 0 ;
    }

  };
   
  // functor that populate compute buffer's extended storage for particle charges
  struct RangeNeighborsClassifier
  {
    double r0min = 0.0;
    double r0max = 0.0;
    double r1min = 0.0;
    double r1max = 0.0;
    double r2min = 0.0;
    double r2max = 0.0;

    template<typename ComputeBufferT, typename FieldArraysT, class NbhDataT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, const FieldArraysT * cells, size_t cell_b, size_t p_b, const NbhDataT& nbh_data) const noexcept
    {
      assert( ssize_t(tab.count) < ssize_t(tab.MaxNeighbors) );
      
      const unsigned int i = tab.count;
      const double r2 = tab.d2[i];
      const double r = sqrt( r2 );
      if( i == 0 )
      {
        tab.ext.nbh_count0 = 0;
        tab.ext.nbh_count1 = 0;
        tab.ext.nbh_count2 = 0;
      }
      unsigned int type = 0;
      unsigned int shift = 0;
      if( i%2 != 0 )
      {
        type = tab.ext.nbh_type[i/2];
        shift = 4;
      }
      
      if( r >= r0min && r < r0max )
      {
        type |= ( 0x01 << shift );
        tab.ext.nbh_list0[ tab.ext.nbh_count0 ++ ] = i;
      }

      if( r >= r1min && r < r1max )
      {
        type |= ( 0x02 << shift );
        tab.ext.nbh_list1[ tab.ext.nbh_count1 ++ ] = i;
      }

      if( r >= r2min && r < r2max )
      {
        type |= ( 0x04 << shift );
        tab.ext.nbh_list2[ tab.ext.nbh_count2 ++ ] = i;
      }
      
      tab.ext.nbh_type[i/2] = type;

      DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, nbh_data );
    }

  };

  struct RangeNeighborsInitComputeBuffer
  {
    const RangeNeighborsClassifier m_classifier;
  
    template<class CPBufT>
    inline void operator () (CPBufT& cpbuf) const
    {
      cpbuf.process_neighbor = m_classifier;
    }
  };

}

