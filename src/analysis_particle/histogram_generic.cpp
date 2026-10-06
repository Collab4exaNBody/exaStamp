/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/histogram.h>
#include <exaStamp/compute/field_combiners.h>

#include <mpi.h>
#include <onika/mpi/data_types.h>

#include <memory>
#include <vector>
#include <string>
#include <type_traits>
#include <limits>

namespace exaStamp
{
  using namespace exanb;

  // HistFieldOrCombiner may either be a plain field tag (field::_ep, ...), wrapped here as a
  // FieldId, or an already-instantiated onika::soatl::FieldCombiner (VelocityNormCombiner, ...),
  // used as-is. Both support cells[i][selector][j] element access.
  template<class T> struct IsFieldCombiner : std::false_type {};
  template<class FuncT, class... fids> struct IsFieldCombiner< onika::soatl::FieldCombiner<FuncT,fids...> > : std::true_type {};

  // combiners whose functor needs the species table at runtime (mass-weighted quantities) can't
  // just be default-constructed like VelocityNormCombiner/ForceNormCombiner can.
  template<class T> struct NeedsSpeciesData : std::false_type {};
  template<> struct NeedsSpeciesData<MonomatKineticEnergyCombiner>  : std::true_type {};
  template<> struct NeedsSpeciesData<MultimatKineticEnergyCombiner> : std::true_type {};
  template<> struct NeedsSpeciesData<MonomatMassCombiner>           : std::true_type {};
  template<> struct NeedsSpeciesData<MultimatMassCombiner>          : std::true_type {};

  template<class HistFieldOrCombiner>
  static inline auto make_hist_field_selector( const ParticleSpecies& species )
  {
    if constexpr ( NeedsSpeciesData<HistFieldOrCombiner>::value )
    {
      return HistFieldOrCombiner{ { species.data() , 0 } };
    }
    else if constexpr ( IsFieldCombiner<HistFieldOrCombiner>::value )
    {
      return HistFieldOrCombiner{};
    }
    else
    {
      return onika::soatl::FieldId<HistFieldOrCombiner>{};
    }
  }

  // single operator, registered once, dispatching at runtime on the "field" slot instead of
  // having one registered operator name per histogrammed quantity
  template<class GridT>
  struct HistogramGenericOperator : public OperatorNode
  {
    // mass-weighted combiners need field::_type to distinguish species when the grid carries
    // one (multi-species run) and fall back to a single-species functor otherwise.
    static constexpr bool has_field_type = GridHasField<GridT,field::_type>::value;
    using KineticEnergyCombinerT = std::conditional_t< has_field_type , MultimatKineticEnergyCombiner , MonomatKineticEnergyCombiner >;
    using MassCombinerT          = std::conditional_t< has_field_type , MultimatMassCombiner           , MonomatMassCombiner >;

    ADD_SLOT( MPI_Comm       , mpi       , INPUT , REQUIRED );
    ADD_SLOT( GridT          , grid      , INPUT , REQUIRED );
    ADD_SLOT( ParticleSpecies, species   , INPUT , REQUIRED );
    ADD_SLOT( std::string    , field     , INPUT , REQUIRED , DocString{"quantity to histogram: ep, charge, vx, vy, vz, fx, fy, fz, rx, ry, rz, vnorm (|v|), fnorm (|f|), vnorm2, fnorm2, mv2 (kinetic energy), mass"} );
    ADD_SLOT( long           , samples   , INPUT , 1000 );
    ADD_SLOT( bool           , ghost     , INPUT , false );
    ADD_SLOT( double         , hist_min  , INPUT , OPTIONAL , DocString{"if set together with hist_max, clamps the histogram interval instead of computing it from the data"} );
    ADD_SLOT( double         , hist_max  , INPUT , OPTIONAL , DocString{"if set together with hist_min, clamps the histogram interval instead of computing it from the data"} );
    ADD_SLOT( Histogram<>    , histogram , OUTPUT );

    inline void execute () override final
    {
      const std::string& f = *field;
           if( f == "ep"       ) { run_histogram<field::_ep>(); }
      else if( f == "charge"   ) { run_histogram<field::_charge>(); }
      else if( f == "vx"       ) { run_histogram<field::_vx>(); }
      else if( f == "vy"       ) { run_histogram<field::_vy>(); }
      else if( f == "vz"       ) { run_histogram<field::_vz>(); }
      else if( f == "fx"       ) { run_histogram<field::_fx>(); }
      else if( f == "fy"       ) { run_histogram<field::_fy>(); }
      else if( f == "fz"       ) { run_histogram<field::_fz>(); }
      else if( f == "rx"       ) { run_histogram<field::_rx>(); }
      else if( f == "ry"       ) { run_histogram<field::_ry>(); }
      else if( f == "rz"       ) { run_histogram<field::_rz>(); }
      else if( f == "vnorm"    ) { run_histogram<VelocityNormCombiner>(); }
      else if( f == "fnorm"    ) { run_histogram<ForceNormCombiner>(); }
      else if( f == "vnorm2"   ) { run_histogram<VelocityNorm2Combiner>(); }
      else if( f == "fnorm2"   ) { run_histogram<ForceNorm2Combiner>(); }
      else if( f == "mv2"      ) { run_histogram<KineticEnergyCombinerT>(); }
      else if( f == "mass"     ) { run_histogram<MassCombinerT>(); }
      else
      {
        lerr << "histogram_generic: unknown field '"<<f<<"'" << std::endl;
        std::abort();
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(
Histograms a single per-particle quantity, selected at runtime via the "field" slot
(one operator instead of one histogram_xxx per quantity).

Usage example:
  - histogram_generic: { field: "ep", samples: 50 }
  - print_histogram: { message: "energy" }
  - write_histogram: { filename: "energy_histogram.csv" }

Supported field values: ep, charge, vx, vy, vz, fx, fy, fz, rx, ry, rz,
vnorm (|v|), fnorm (|f|), vnorm2 (|v|^2), fnorm2 (|f|^2), mv2 (kinetic energy), mass.
mv2/mass automatically use the multi- or mono-species functor depending on whether
this grid variant carries a per-particle "type" field.
hist_min/hist_max optionally clamp the interval instead of computing it from the data.
)EOF";
    }

  private:

    template<class HistFieldOrCombiner>
    inline void run_histogram()
    {
      if constexpr ( ! GridHasField<GridT,HistFieldOrCombiner>::value )
      {
        lerr << "histogram_generic: field '"<<(*field)<<"' is not available on this grid" << std::endl;
        std::abort();
      }
      else
      {
        using FieldSelector = decltype( make_hist_field_selector<HistFieldOrCombiner>(*species) );
        using ValueType = typename FieldSelector::value_type;
        const FieldSelector hist_field = make_hist_field_selector<HistFieldOrCombiner>(*species);

        MPI_Comm comm = *mpi;
        int nprocs = 1;
        int rank = 0;
        MPI_Comm_size(comm,&nprocs);
        MPI_Comm_rank(comm,&rank);

        auto cells = grid->cells();
        IJK dims = grid->dimension();
        ssize_t gl = grid->ghost_layers();
        if( *ghost ) { gl = 0; }

        // min max computation, unless a clamp interval was provided
        ValueType min_val;
        ValueType max_val;

        if( hist_min.has_value() && hist_max.has_value() )
        {
          min_val = *hist_min;
          max_val = *hist_max;
        }
        else
        {
          min_val = std::numeric_limits<ValueType>::max();
          max_val = std::numeric_limits<ValueType>::lowest();

#         pragma omp parallel
          {
            ValueType local_min_val = std::numeric_limits<ValueType>::max();
            ValueType local_max_val = std::numeric_limits<ValueType>::lowest();

            GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc)
            {
              size_t i = grid_ijk_to_index( dims , loc + gl );
              size_t n = cells[i].size();
              auto value_ptr = cells[i][hist_field];

              for(size_t j=0;j<n;j++)
              {
                ValueType x = value_ptr[j];
                local_min_val = std::min( local_min_val , x );
                local_max_val = std::max( local_max_val , x );
              }
            }
            GRID_OMP_FOR_END
#           pragma omp critical
            {
              min_val = std::min( local_min_val , min_val );
              max_val = std::max( local_max_val , max_val );
            }
          }

          // MPI min/max
          if( nprocs > 1 )
          {
            ValueType tmp[2] = { -min_val , max_val };
            MPI_Allreduce(MPI_IN_PLACE,tmp,2, onika::mpi::mpi_datatype<ValueType>() ,MPI_MAX,comm);
            min_val = - tmp[0];
            max_val = tmp[1];
          }

          // a single bound may be provided while the other is auto-computed
          if( hist_min.has_value() ) { min_val = *hist_min; }
          if( hist_max.has_value() ) { max_val = *hist_max; }
        }

        // histogram counting
        size_t hist_size = *samples;
        histogram->m_min_val = min_val;
        histogram->m_max_val = max_val;
        histogram->m_data.resize( hist_size );

        int max_nt = omp_get_max_threads();
        double* per_thread_histogram[max_nt];

#       pragma omp parallel
        {
          int nt = omp_get_num_threads();
          int tid = omp_get_thread_num();
          assert( tid<max_nt && nt<=max_nt );

          // stack allocated, per thread histogram
          double local_hist[hist_size];
          for(size_t i=0;i<hist_size;i++) { local_hist[i] = 0.0; }

          GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc)
          {
            size_t i = grid_ijk_to_index( dims , loc + gl );
            auto value_ptr = cells[i][hist_field];
            size_t n = cells[i].size();
            for(size_t j=0;j<n;j++)
            {
              ValueType v = value_ptr[j];
              ssize_t bin = static_cast<size_t>( ( (v-min_val) * hist_size ) / ( max_val - min_val ) );
              if( bin < 0 ) { bin=0; }
              if( bin >= static_cast<ssize_t>(hist_size) ) { bin = hist_size-1; }
              local_hist[bin] += 1.0;
            }
          }
          GRID_OMP_FOR_END

          per_thread_histogram[tid] = local_hist;
          size_t start = ( hist_size * tid ) / nt;
          size_t end = ( hist_size * (tid+1) ) / nt;

#         pragma omp barrier

          double* h = per_thread_histogram[0];
          for(size_t i=start;i<end;i++)
          {
            histogram->m_data[i] = h[i];
          }
          for(int t=1;t<nt;t++)
          {
            h = per_thread_histogram[t];
            for(size_t i=start;i<end;i++)
            {
              histogram->m_data[i] += h[i];
            }
          }
        }

        if( nprocs > 1 )
        {
          MPI_Allreduce(MPI_IN_PLACE,histogram->m_data.data(),hist_size,MPI_DOUBLE,MPI_SUM,comm);
        }
      }
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(histogram_generic)
  {
    OperatorNodeFactory::instance()->register_factory( "histogram_generic" , make_grid_variant_operator< HistogramGenericOperator > );
  }

}
