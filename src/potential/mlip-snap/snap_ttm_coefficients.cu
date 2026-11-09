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

// Electronic temperature dependent SNAP coefficients, equivalent of LAMMPS pair_style snap_ttm
// (ML-SNAP/pair_snap_ttm.cpp). For each particle, Te is read from the grid_cell_values "te" field
// (subcell containing the particle, as LAMMPS fix ttm does) and converted to eV, then :
//   beta0(Te)            = polynomial (bzeropolyorder), Horner, highest order first
//   beta1..ncoeff(Te)    = natural cubic spline through the tabulated values (rows 1..ncoeffall-1 of
//                          the betas file, row 0 is unused), cubic extrapolation of the end intervals
// Result is written in snap_atom_coefs ( ncoeffall values per particle, layout of one .snapcoeff element
// block, indexed by cell_particle_offset[cell]+p ), consumed by snap_force_fp64.
//
// betas file format (same as LAMMPS) :
//   ntelec N
//   Te_0 ... Te_N-1                 (eV, strictly increasing)
//   bzeropolyorder P
//   c_P ... c_0                     (P+1 values, highest order first)
//   ncoeffall lines of N values     (beta_r at each Te knot)

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/memory/allocator.h>
#include <onika/file_utils.h>
#include <onika/log.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/grid_cell_particles/grid_cell_values_utils.h>

#include <mpi.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

namespace exaStamp
{
  using namespace exanb;

  // spline table : knots, per row and interval cubic coefficients (a,b,c,d), beta0 polynomial
  struct SnapTtmBetaTable
  {
    int N = 0;   // number of Te knots
    int M = 0;   // number of rows (ncoeffall)
    int P = -1;  // beta0 polynomial order
    onika::memory::CudaMMVector<double> knots;   // N
    onika::memory::CudaMMVector<double> coeffs;  // M * (N-1) * 4
    onika::memory::CudaMMVector<double> poly;    // P+1
    std::vector< std::vector<double> > betas;    // M x N, raw values (host only)
  };

  struct SnapTtmCoefsFunctor
  {
    const size_t * __restrict__ cell_particle_offset = nullptr;
    double * __restrict__ out = nullptr;

    const double * __restrict__ knots = nullptr;
    const double * __restrict__ coeffs = nullptr;
    const double * __restrict__ poly = nullptr;
    int N = 0;
    int M = 0;
    int P = -1;

    const double * __restrict__ te_ptr = nullptr;
    size_t te_stride = 0;
    IJK grid_dims = { 0, 0, 0 };
    IJK grid_offset = { 0, 0, 0 };
    Vec3d grid_origin = { 0., 0., 0. };
    double cell_size = 0.0;
    double subcell_size = 0.0;
    ssize_t subdiv = 1;
    double te_conv = 0.0; // Te (K) -> Te (eV)

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t cell_i, unsigned int p_i, double rx, double ry, double rz ) const
    {
      using namespace GridCellValuesUtils;

      // Te of the subcell containing the particle (LAMMPS fix ttm: nearest grid cell, no interpolation).
      // a particle may sit slightly outside its own cell between two grid updates, hence the neighbor cell.
      const IJK cell_loc = grid_index_to_ijk( grid_dims, cell_i );
      const Vec3d cell_origin = grid_origin + ( (grid_offset + cell_loc) * cell_size );
      const Vec3d rco = Vec3d{ rx, ry, rz } - cell_origin;
      IJK te_cell_loc, te_subcell_loc;
      localize_subcell( rco, cell_size, subcell_size, subdiv, te_cell_loc, te_subcell_loc );
      te_cell_loc = te_cell_loc + cell_loc;
      te_cell_loc.i = onika::cuda::max( ssize_t(0) , onika::cuda::min( te_cell_loc.i , grid_dims.i - 1 ) );
      te_cell_loc.j = onika::cuda::max( ssize_t(0) , onika::cuda::min( te_cell_loc.j , grid_dims.j - 1 ) );
      te_cell_loc.k = onika::cuda::max( ssize_t(0) , onika::cuda::min( te_cell_loc.k , grid_dims.k - 1 ) );
      const ssize_t te_cell_i = grid_ijk_to_index( grid_dims, te_cell_loc );
      const ssize_t te_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv}, te_subcell_loc );
      const double te = te_ptr[ te_cell_i * te_stride + te_subcell_i ] * te_conv;

      double * __restrict__ coefi = out + size_t(M) * ( cell_particle_offset[cell_i] + p_i );

      // beta0 : polynomial, Horner
      double b0 = 0.0;
      for(int i=0;i<=P;i++) b0 = b0 * te + poly[i];
      coefi[0] = b0;

      // interval k such that knots[k] <= te < knots[k+1], clamped to [0,N-2] (extrapolation outside)
      int k = 0;
      while( k < N-2 && te >= knots[k+1] ) ++k;
      const double t = te - knots[k];
      for(int row=1;row<M;row++)
      {
        const double * __restrict__ c = coeffs + ( size_t(row) * (N-1) + k ) * 4;
        coefi[row] = ( ( c[3] * t + c[2] ) * t + c[1] ) * t + c[0];
      }
    }
  };
}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits<exaStamp::SnapTtmCoefsFunctor>
  {
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  // natural cubic spline, straight port of LAMMPS PairSNAPTTM::factor_tridiagonal_natural,
  // solve_tridiagonal_natural and build_row_coeffs (bit-for-bit same operation order)
  static inline void snap_ttm_factor_tridiagonal_natural( const std::vector<double>& h, std::vector<double>& cprime, std::vector<double>& denom )
  {
    const int N = (int)denom.size();
    std::vector<double> a(N, 0.0), b(N, 0.0), c(N, 0.0);
    b[0] = 1.0;
    b[N-1] = 1.0;
    for (int i = 1; i <= N-2; ++i) { a[i] = h[i-1]; b[i] = 2.0 * (h[i-1] + h[i]); c[i] = h[i]; }
    cprime[0] = (N > 1) ? c[0] / b[0] : 0.0;
    denom[0] = b[0];
    for (int i = 1; i < N; ++i)
    {
      denom[i] = b[i] - a[i] * cprime[i-1];
      cprime[i] = (i == N-1) ? 0.0 : c[i] / denom[i];
    }
  }

  static inline void snap_ttm_solve_tridiagonal_natural( const std::vector<double>& h, const std::vector<double>& cprime, const std::vector<double>& denom, std::vector<double>& rhs )
  {
    const int N = (int)rhs.size();
    std::vector<double> y(N);
    y[0] = rhs[0] / denom[0];
    for (int i = 1; i < N; ++i) y[i] = (rhs[i] - h[i-1] * y[i-1]) / denom[i];
    rhs[N-1] = y[N-1];
    for (int i = N - 2; i >= 0; --i) rhs[i] = y[i] - cprime[i] * rhs[i+1];
  }

  static inline void snap_ttm_build_table( SnapTtmBetaTable& tab )
  {
    const int N = tab.N;
    const int M = tab.M;
    std::vector<double> h( N-1 );
    for(int i=1;i<N;i++) h[i-1] = tab.knots[i] - tab.knots[i-1];
    std::vector<double> cprime( N, 0.0 ), denom( N, 0.0 );
    snap_ttm_factor_tridiagonal_natural( h, cprime, denom );
    tab.coeffs.assign( size_t(M) * (N-1) * 4 , 0.0 );
    for(int row=0;row<M;row++)
    {
      const double * y = tab.betas[row].data();
      std::vector<double> rhs( N, 0.0 );
      for (int i = 1; i <= N-2; ++i)
      {
        const double slope_next = (y[i+1] - y[i]) / h[i];
        const double slope_prev = (y[i] - y[i-1]) / h[i-1];
        rhs[i] = 3.0 * (slope_next - slope_prev);
      }
      snap_ttm_solve_tridiagonal_natural( h, cprime, denom, rhs );
      const std::vector<double>& cvec = rhs;
      for (int k = 0; k < N - 1; ++k)
      {
        const double ak = y[k];
        const double ck = cvec[k];
        const double dk = (cvec[k+1] - cvec[k]) / (3.0 * h[k]);
        const double bk = (y[k+1] - y[k]) / h[k] - (2.0*ck + cvec[k+1]) * h[k] / 3.0;
        const size_t base = ( size_t(row) * (N-1) + k ) * 4;
        tab.coeffs[base+0] = ak;
        tab.coeffs[base+1] = bk;
        tab.coeffs[base+2] = ck;
        tab.coeffs[base+3] = dk;
      }
    }
  }

  template< class GridT >
  class SnapTtmCoefficients : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi              , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT           , grid             , INPUT , REQUIRED );
    ADD_SLOT( Domain          , domain           , INPUT , REQUIRED );
    ADD_SLOT( GridCellValues  , grid_cell_values , INPUT , OPTIONAL , DocString{"must hold the electronic temperature field 'te' (K). optional only because preinit_rcut_max runs compute_force on an empty grid"} );
    ADD_SLOT( std::string     , betas_file       , INPUT , REQUIRED , DocString{"LAMMPS snap_ttm betas file (.snapbetas)"} );
    ADD_SLOT( double          , te_conv_factor   , INPUT , 8.61732814974056e-5 , DocString{"Te (K) to betas file Te unit (eV) factor. default is LAMMPS pair snap_ttm's hard coded value"} );
    ADD_SLOT( std::string     , betas_check_file , INPUT , OPTIONAL , DocString{"if set, writes beta0..betaM-1 evaluated on 10000 Te points in [0,6] eV to this file (rank 0, first call)"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , snap_atom_coefs , INPUT_OUTPUT , DocString{"per-particle SNAP coefficients, input of snap_force_fp64"} );
    ADD_SLOT( SnapTtmBetaTable , snap_ttm_table  , PRIVATE );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Computes per-particle electronic temperature dependent SNAP coefficients (LAMMPS pair_style snap_ttm),
from the grid_cell_values 'te' field. To be placed right before snap_force_fp64 in compute_force.
)EOF";
    }

    inline void execute () override final
    {
      if( snap_ttm_table->N == 0 ) read_betas_file();

      auto & tab = *snap_ttm_table;
      const size_t total_particles = grid->number_of_particles();
      snap_atom_coefs->resize( total_particles * tab.M );
      if( grid->number_of_cells() == 0 ) return;

      if( ! grid_cell_values.has_value() || ! grid_cell_values->has_field("te") )
      {
        fatal_error() << pathname() << ": grid_cell_values has no 'te' field (use init_ttm in setup_system)" << std::endl;
      }
      const ssize_t subdiv = grid_cell_values->field("te").m_subdiv;
      auto cell_te_data = grid_cell_values->field_data("te");
      const double cell_size = domain->cell_size();

      SnapTtmCoefsFunctor func = {
        grid->cell_particle_offset_data(), snap_atom_coefs->data(),
        tab.knots.data(), tab.coeffs.data(), tab.poly.data(), tab.N, tab.M, tab.P,
        cell_te_data.m_data_ptr, cell_te_data.m_stride,
        grid->dimension(), grid->offset(), grid->origin(),
        cell_size, cell_size / subdiv, subdiv,
        *te_conv_factor };
      compute_cell_particles( *grid, true, func, FieldSet<field::_rx,field::_ry,field::_rz>{}, parallel_execution_context() );
    }

  private:
    inline void read_betas_file()
    {
      auto & tab = *snap_ttm_table;
      int rank = 0;
      MPI_Comm_rank( *mpi, &rank );

      // rank 0 reads, flat buffer [ N, P, M, knots(N), poly(P+1), betas(M*N) ] broadcast to others
      std::vector<double> buf;
      if( rank == 0 )
      {
        const std::string file_name = onika::data_file_path( *betas_file );
        std::ifstream fin( file_name );
        if( ! fin ) fatal_error() << pathname() << ": cannot open betas file " << file_name << std::endl;
        std::vector< std::vector<std::string> > lines;
        std::string line;
        while( std::getline(fin,line) )
        {
          line = line.substr( 0, line.find('#') );
          std::istringstream iss(line);
          std::vector<std::string> words; std::string w;
          while( iss >> w ) words.push_back(w);
          if( ! words.empty() ) lines.push_back( words );
        }
        if( lines.size() < 5 || lines[0].size() != 2 || lines[2].size() != 2 )
          fatal_error() << pathname() << ": bad betas file format in " << file_name << std::endl;
        const int N = std::stoi( lines[0][1] );
        const int P = std::stoi( lines[2][1] );
        const int M = lines.size() - 4;
        if( N < 2 || int(lines[1].size()) != N || P < 0 || int(lines[3].size()) != P+1 )
          fatal_error() << pathname() << ": bad betas file format in " << file_name << std::endl;
        buf.push_back(N); buf.push_back(P); buf.push_back(M);
        for(int i=0;i<N;i++) buf.push_back( std::stod(lines[1][i]) );
        for(int i=0;i<=P;i++) buf.push_back( std::stod(lines[3][i]) );
        for(int r=0;r<M;r++)
        {
          if( int(lines[4+r].size()) != N ) fatal_error() << pathname() << ": betas row "<<r<<" has "<<lines[4+r].size()<<" values, expected "<<N<< std::endl;
          for(int i=0;i<N;i++) buf.push_back( std::stod(lines[4+r][i]) );
        }
      }
      long bufsize = buf.size();
      MPI_Bcast( &bufsize, 1, MPI_LONG, 0, *mpi );
      buf.resize( bufsize );
      MPI_Bcast( buf.data(), bufsize, MPI_DOUBLE, 0, *mpi );

      size_t pos = 0;
      tab.N = buf[pos++]; tab.P = buf[pos++]; tab.M = buf[pos++];
      tab.knots.assign( buf.begin()+pos, buf.begin()+pos+tab.N ); pos += tab.N;
      tab.poly.assign( buf.begin()+pos, buf.begin()+pos+tab.P+1 ); pos += tab.P+1;
      tab.betas.assign( tab.M, std::vector<double>(tab.N) );
      for(int r=0;r<tab.M;r++) for(int i=0;i<tab.N;i++) tab.betas[r][i] = buf[pos++];
      for(int i=1;i<tab.N;i++)
      {
        if( !( tab.knots[i] > tab.knots[i-1] ) ) fatal_error() << pathname() << ": betas file Te knots must be strictly increasing" << std::endl;
      }

      snap_ttm_build_table( tab );
      ldbg << pathname() << ": " << tab.N << " Te knots ["<<tab.knots[0]<<","<<tab.knots[tab.N-1]<<"] eV, "<< tab.M << " coefficients, beta0 polynomial order " << tab.P << std::endl;

      if( rank == 0 && betas_check_file.has_value() ) write_check_file();
    }

    inline void write_check_file()
    {
      const auto & tab = *snap_ttm_table;
      std::ofstream fout( *betas_check_file );
      fout << "# Te b0(Te) b1(Te) ... b" << (tab.M-1) << "(Te)\n" << std::fixed << std::setprecision(17);
      const int npts = 10000;
      std::vector<double> coefi( tab.M );
      for(int q=0;q<npts;q++)
      {
        const double te = 6.0 * static_cast<double>(q) / static_cast<double>(npts - 1);
        double b0 = 0.0;
        for(int i=0;i<=tab.P;i++) b0 = b0 * te + tab.poly[i];
        int k = 0;
        while( k < tab.N-2 && te >= tab.knots[k+1] ) ++k;
        const double t = te - tab.knots[k];
        fout << te << " " << b0;
        for(int row=1;row<tab.M;row++)
        {
          const double * c = tab.coeffs.data() + ( size_t(row) * (tab.N-1) + k ) * 4;
          fout << " " << ( ( c[3] * t + c[2] ) * t + c[1] ) * t + c[0];
        }
        fout << "\n";
      }
    }
  };

  template<class GridT> using SnapTtmCoefficientsTmpl = SnapTtmCoefficients<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(snap_ttm_coefficients)
  {
    OperatorNodeFactory::instance()->register_factory( "snap_ttm_coefficients" , make_grid_variant_operator< SnapTtmCoefficientsTmpl > );
  }

}
