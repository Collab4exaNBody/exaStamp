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

#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/mlip_ye3t/ye3t_context.h>
#include <exaStamp/mlip_ye3t/ye3t_math.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/file_utils.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace exaStamp
{
  using namespace exanb;
  using namespace mlip_ye3t;

  class Ye3tInit : public OperatorNode
  {
    ADD_SLOT( Ye3tParams      , parameters , INPUT        , REQUIRED , DocString{"file: .yace model (linear ACE, ye3t-lammps pair_style ye3t), or .json tagged-Cauchy composite / slice v2; analytic, mem_gb, emulate: see documentation"} );
    ADD_SLOT( double          , rcut_max   , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( ParticleSpecies , species    , INPUT        , REQUIRED );
    ADD_SLOT( Ye3tContext     , ye3t_ctx   , OUTPUT );

    std::string m_fname;

    template<class... T> [[noreturn]] void fail(const T&... msg)
    {
      std::ostringstream oss;
      ( oss << ... << msg );
      fatal_error() << "ye3t_init: " << oss.str() << " (" << m_fname << ")" << std::endl;
      std::abort();
    }

    YAML::Node req(const YAML::Node& n, const char* key, const std::string& where)
    {
      const YAML::Node x = n[key];
      if( !x ) fail(where, ": missing '", key, "'");
      return x;
    }

    // spline coefficients [(n+1)][nf][4] of nf functions with values v / derivatives d at the nodes r = h (k+1), k < n
    static void build_spline(const std::vector<double>& v, const std::vector<double>& d, int n, int nf, double h, double* c)
    {
      std::fill( c, c + size_t(n + 1) * nf * 4, 0.0 );
      for(int k=1;k<=n;k++)
      {
        const bool right = k < n;
        for(int f=0;f<nf;f++)
        {
          const double f0 = v[size_t(k-1)*nf+f], f1 = right ? v[size_t(k)*nf+f] : 0.0;
          const double d0 = d[size_t(k-1)*nf+f] * h, d1 = right ? d[size_t(k)*nf+f] * h : 0.0;
          double* o = c + ( size_t(k) * nf + f ) * 4;
          o[0] = f0;
          o[1] = d0;
          o[2] = 3.0 * ( f1 - f0 ) - d1 - 2.0 * d0;
          o[3] = -2.0 * ( f1 - f0 ) + d1 + d0;
        }
      }
    }

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(

Linear atomic cluster expansion from a PACE-compatible .yace file (the .yace family of ye3t-lammps pair_style ye3t:
ChebExpCos radial basis, one linear density, no hard core, no inner cutoff), for ye3t_force. Native onika kernels
(GPU if available, else OpenMP).

ye3t_init:
  parameters:
    file: model.yace
    analytic: false   # true: radial functions evaluated directly instead of the cubic splines of the model (debug)
    mem_gb: 1.0       # budget of the per-chunk buffers
    emulate: 0        # debug: N > 0 runs the GPU code paths on the host with blocks of N threads (slow)

Tagged-Cauchy models (ye3t-lammps model_family tagged_cauchy, slice v2): file is a .json, either a composite
manifest (schema ye3t_tagged_cauchy_composite_v1: a .yace backbone plus a tagged correction, paths relative to the
manifest, e.g. model.ye3t.json of the ye3t-lammps cost_comparison models) or a tagged slice alone (schema
ye3t_tagged_cauchy_slice_v2). The correction adds to E_i densities and same-edge moments of real tesseral components
times a falling factorial of the edge count; its radial functions are always evaluated directly. Models fitted as a
residual to ZBL need the ZBL pair term of the fit (pair_style hybrid/overlay ye3t ... zbl inner outer in LAMMPS):

compute_force:
  - ye3t_force
  - zbl_multi_force:
      rcut: 1.649 ang
      parameters:
        - { type_a: Cu , type_b: Cu , rcut: 1.649 ang , parameters: { r1: 1.283 ang , rc: 1.649 ang } }

Physical-image tagged models (slices v3 / v4) and lifted-Cauchy models are not supported.

The exaStamp species are mapped to the model elements by name.

)EOF";
    }

    // .yace model (also the backbone of a tagged-Cauchy composite)
    std::shared_ptr<Ye3tTables> load_yace( const std::string& path )
    {
      m_fname = path;
      YAML::Node root;
      try { root = YAML::LoadFile( m_fname ); }
      catch( const YAML::Exception& e ) { fail("cannot read model file: ", e.what()); }

      auto T = std::make_shared<Ye3tTables>();
      T->analytic = parameters->analytic;

      // elements, energies
      const auto elements = req(root, "elements", "model").as< std::vector<std::string> >();
      const int nel = elements.size();
      if( nel == 0 ) fail("no element");
      T->elements = elements;
      T->nel = nel;
      const auto e0 = req(root, "E0", "model").as< std::vector<double> >();
      if( int(e0.size()) != nel ) fail("E0 has ", e0.size(), " values for ", nel, " elements");
      const double spacing = req(root, "deltaSplineBins", "model").as<double>();
      if( !( spacing > 0.0 ) ) fail("deltaSplineBins must be positive");
      T->e0.assign( e0.begin(), e0.end() );
      T->fs0.assign( nel, 0.0 );
      T->rho_lim.assign( nel, 0.0 );
      const YAML::Node emb = req(root, "embeddings", "model");
      for(int mu=0;mu<nel;mu++)
      {
        const YAML::Node e = emb[mu];
        const std::string w = "embeddings." + std::to_string(mu);
        if( !e ) fail(w, ": missing");
        if( req(e, "ndensity", w).as<int>() != 1 ) fail(w, ": only one density is supported");
        const auto fs = req(e, "FS_parameters", w).as< std::vector<double> >();
        if( fs.size() != 2 || fs[1] != 1.0 ) fail(w, ": only the linear embedding (FS_parameters [c, 1]) is supported");
        const std::string npoti = req(e, "npoti", w).as<std::string>();
        if( npoti != "FinnisSinclair" && npoti != "FinnisSinclairShiftedScaled" ) fail(w, ": unsupported embedding ", npoti);
        const double rc = req(e, "rho_core_cutoff", w).as<double>(), drc = req(e, "drho_core_cutoff", w).as<double>();
        if( !( rc > 0.0 ) || drc < 0.0 || drc >= rc ) fail(w, ": invalid density core cutoff");
        T->fs0[mu] = fs[0];
        T->rho_lim[mu] = rc - drc;
      }

      // bonds, radial splines
      const YAML::Node bonds = req(root, "bonds", "model");
      if( !bonds.IsMap() || int(bonds.size()) != nel * nel ) fail("bonds: expected every ordered element pair");
      T->bonds.resize( nel * nel );
      std::vector<double> pool;
      std::vector<bool> seen( nel * nel, false );
      for(const auto& kv : bonds)
      {
        const auto key = kv.first.as< std::vector<int> >();
        if( key.size() != 2 || key[0] < 0 || key[0] >= nel || key[1] < 0 || key[1] >= nel ) fail("bonds: bad key");
        const int bi = key[0] * nel + key[1];
        if( seen[bi] ) fail("bonds: duplicate pair");
        seen[bi] = true;
        const YAML::Node n = kv.second;
        const std::string w = "bonds [" + std::to_string(key[0]) + "," + std::to_string(key[1]) + "]";
        if( req(n, "radbasename", w).as<std::string>() != "ChebExpCos" ) fail(w, ": only ChebExpCos is supported");
        if( req(n, "prehc", w).as<double>() != 0.0 ) fail(w, ": hard-core repulsion is not supported");
        if( req(n, "rcut_in", w).as<double>() != 0.0 || req(n, "dcut_in", w).as<double>() != 0.0 ) fail(w, ": inner cutoffs are not supported");
        if( n["inner_cutoff_type"] && n["inner_cutoff_type"].as<std::string>() != "distance" ) fail(w, ": unsupported inner_cutoff_type");
        Ye3tBond b;
        b.nrad = req(n, "nradmax", w).as<int>();
        b.lb = req(n, "lmax", w).as<int>();
        b.nbase = req(n, "nradbasemax", w).as<int>();
        const auto rp = req(n, "radparameters", w).as< std::vector<double> >();
        if( rp.empty() || !( rp[0] > 0.0 ) ) fail(w, ": radparameters[0] must be positive");
        b.lambda = rp[0];
        b.rc = req(n, "rcut", w).as<double>();
        b.dcut = req(n, "dcut", w).as<double>();
        if( !( b.rc > 0.0 ) || !( b.dcut > 0.0 ) || b.dcut >= b.rc ) fail(w, ": invalid outer cutoff");
        if( b.nrad < 1 || b.nbase < 1 || b.lb < 0 ) fail(w, ": invalid sizes");
        if( b.nbase > 64 ) fail(w, ": at most 64 radial base functions");
        if( b.lb > kYe3tMaxLmax ) fail(w, ": lmax ", b.lb, " > ", kYe3tMaxLmax);
        b.ncontr = b.nrad * ( b.lb + 1 );
        const auto rc3 = req(n, "radcoefficients", w).as< std::vector< std::vector< std::vector<double> > > >();
        if( int(rc3.size()) != b.nrad ) fail(w, ": radcoefficients: wrong radial dimension");
        b.crad = pool.size();
        for(int nn=0;nn<b.nrad;nn++)
        {
          if( int(rc3[nn].size()) != b.lb + 1 ) fail(w, ": radcoefficients: wrong angular dimension");
          for(int l=0;l<=b.lb;l++)
          {
            if( int(rc3[nn][l].size()) != b.nbase ) fail(w, ": radcoefficients: wrong base dimension");
            pool.insert( pool.end(), rc3[nn][l].begin(), rc3[nn][l].end() );
          }
        }
        // spline tables
        const double ratio = b.rc / spacing;
        if( ratio < 1.0 || ratio > 1e8 ) fail(w, ": spline interval count out of range");
        b.nint = int( ratio );
        b.scale = double( b.nint ) / b.rc;
        const double h = b.rc / double( b.nint );
        std::vector<double> gv( size_t(b.nint) * b.nbase ), gd( gv.size() );
        for(int k=0;k<b.nint;k++)
        {
          const double r = h * double( k + 1 );
          double* g = gv.data() + size_t(k) * b.nbase;
          double* dg = gd.data() + size_t(k) * b.nbase;
          if( r >= b.rc ) { std::fill( g, g + b.nbase, 0.0 ); std::fill( dg, dg + b.nbase, 0.0 ); }
          else cheb_exp_cos( b, r, g, dg );
        }
        b.bspl = pool.size();
        pool.resize( pool.size() + size_t(b.nint + 1) * b.nbase * 4 );
        build_spline( gv, gd, b.nint, b.nbase, h, pool.data() + b.bspl );
        T->bonds[bi] = b;
        T->lmax = std::max( T->lmax , b.lb );
        T->nbase_max = std::max( T->nbase_max , b.nbase );
        T->ncontr_max = std::max( T->ncontr_max , b.ncontr );
        T->rcut = std::max( T->rcut , b.rc );
      }
      T->pool.assign( pool.begin(), pool.end() );
      T->nlh = ( T->lmax + 1 ) * ( T->lmax + 2 ) / 2;
      T->ylm_a.resize( T->nlh );
      T->ylm_b.resize( T->nlh );
      ylm_coefficients( T->lmax, T->ylm_a.data(), T->ylm_b.data() );

      // functions -> sources, channels, merged monomials, per central element
      const YAML::Node funcs = req(root, "functions", "model");
      std::vector<int> src_kind, src_r, src_lm, src_nu, ch_src, ch_tr, node_ch, node_depth;
      std::vector<double> node_coef;
      std::vector<int> src0{ 0 }, ch0{ 0 }, no0{ 0 }, srcb;
      size_t nmono = 0, nnode = 0;
      size_t nfunc = 0, nrows = 0;
      for(int mu=0;mu<nel;mu++)
      {
        const YAML::Node fl = funcs[mu];
        const std::string wb = "functions." + std::to_string(mu);
        if( !fl || !fl.IsSequence() || fl.size() == 0 ) fail(wb, ": expected a list of functions");
        using CKey = std::tuple<int,int,int,int,int>;    // kind, nu, n, l, m
        std::map<CKey,int> chan;
        std::vector<CKey> chan_keys;
        std::map< std::vector<std::pair<int,int>> , int > mono;
        std::vector< std::vector<std::pair<int,int>> > mono_keys;
        std::vector<double> mono_coef;
        for(size_t fi=0;fi<fl.size();fi++)
        {
          const YAML::Node f = fl[fi];
          const std::string w = wb + "[" + std::to_string(fi) + "]";
          if( req(f, "mu0", w).as<int>() != mu ) fail(w, ": mu0 does not match the element");
          if( req(f, "ndensity", w).as<int>() != 1 ) fail(w, ": only one density is supported");
          const int rank = req(f, "rank", w).as<int>();
          const int nrow = req(f, "num_ms_combs", w).as<int>();
          const auto mus = req(f, "mus", w).as< std::vector<int> >();
          const auto ns = req(f, "ns", w).as< std::vector<int> >();
          const auto ls = req(f, "ls", w).as< std::vector<int> >();
          const auto ms = req(f, "ms_combs", w).as< std::vector<int> >();
          const auto ct = req(f, "ctildes", w).as< std::vector<double> >();
          if( rank < 1 || int(mus.size()) != rank || int(ns.size()) != rank || int(ls.size()) != rank ) fail(w, ": inconsistent rank");
          if( nrow < 1 || int(ms.size()) != rank * nrow || int(ct.size()) != nrow ) fail(w, ": inconsistent ms_combs / ctildes");
          if( rank == 1 && ( ls[0] != 0 || ms[0] != 0 || nrow != 1 ) ) fail(w, ": rank-one functions need one l = m = 0 row");
          ++ nfunc;
          for(int row=0;row<nrow;row++)
          {
            ++ nrows;
            int msum = 0;
            std::vector<int> rc;
            for(int s=0;s<rank;s++)
            {
              const int nu = mus[s], n = ns[s] - 1, l = ls[s], m = ms[size_t(row)*rank+s];
              if( nu < 0 || nu >= nel ) fail(w, ": neighbour element out of range");
              const Ye3tBond& b = T->bonds[ mu * nel + nu ];
              if( n < 0 || n >= ( rank == 1 ? b.nbase : b.nrad ) ) fail(w, ": radial index out of range");
              if( l < 0 || l > b.lb || std::abs(m) > l ) fail(w, ": angular index out of range");
              msum += m;
              const CKey k{ rank == 1 ? 0 : 1, nu, n, l, m };
              auto it = chan.find( k );
              if( it == chan.end() ) { it = chan.emplace( k, int(chan_keys.size()) ).first; chan_keys.push_back( k ); }
              rc.push_back( it->second );
            }
            if( msum != 0 ) fail(w, ": a ms row does not couple to M = 0");
            std::sort( rc.begin(), rc.end() );
            std::vector<std::pair<int,int>> mk;
            for(size_t p=0;p<rc.size();)
            {
              size_t e = p + 1;
              while( e < rc.size() && rc[e] == rc[p] ) ++ e;
              mk.emplace_back( rc[p], int(e - p) );
              p = e;
            }
            auto it = mono.find( mk );
            if( it == mono.end() ) { it = mono.emplace( mk, int(mono_keys.size()) ).first; mono_keys.push_back( mk ); mono_coef.push_back( 0.0 ); }
            mono_coef[ it->second ] += ct[row];
          }
        }
        // sources (m >= 0), sorted by neighbour element, in order of first use
        std::map<CKey,int> srcmap;
        std::vector<CKey> skeys;
        for(const auto& k : chan_keys)
        {
          const CKey s{ std::get<0>(k), std::get<1>(k), std::get<2>(k), std::get<3>(k), std::abs( std::get<4>(k) ) };
          if( srcmap.find( s ) == srcmap.end() ) { srcmap.emplace( s, 0 ); skeys.push_back( s ); }
        }
        std::stable_sort( skeys.begin(), skeys.end(), [](const CKey& a, const CKey& b) { return std::get<1>(a) < std::get<1>(b); } );
        for(size_t q=0;q<skeys.size();q++) srcmap[ skeys[q] ] = int(q);
        for(int nu=0;nu<=nel;nu++)
        {
          int q = 0;
          while( q < int(skeys.size()) && std::get<1>( skeys[q] ) < nu ) ++ q;
          srcb.push_back( q );
        }
        for(const auto& s : skeys)
        {
          const int kind = std::get<0>(s), nu = std::get<1>(s), n = std::get<2>(s), l = std::get<3>(s), m = std::get<4>(s);
          const Ye3tBond& b = T->bonds[ mu * nel + nu ];
          src_kind.push_back( kind );
          src_nu.push_back( nu );
          src_r.push_back( kind == 0 ? n : n * ( b.lb + 1 ) + l );
          src_lm.push_back( kind == 0 ? 0 : lmh( l, m ) );
        }
        for(const auto& k : chan_keys)
        {
          const int m = std::get<4>(k);
          const CKey s{ std::get<0>(k), std::get<1>(k), std::get<2>(k), std::get<3>(k), std::abs(m) };
          ch_src.push_back( srcmap[s] );
          ch_tr.push_back( m < 0 ? ( ( -m ) % 2 == 0 ? 1 : -1 ) : 0 );
        }
        // complex-conjugate symmetry: the monomial with the partner channels (m -> -m) has the conjugate value times
        // prod (s_c / s_partner)^p (channel = s X, X = A or conj A), the same real part: the coefficients must agree
        // (real density, checked as ye3t-lammps does) and one monomial of each pair is kept with the summed coefficient
        std::vector<int> partner( chan_keys.size(), -1 );
        for(size_t c=0;c<chan_keys.size();c++)
        {
          const auto& k = chan_keys[c];
          const auto it = chan.find( CKey{ std::get<0>(k), std::get<1>(k), std::get<2>(k), std::get<3>(k), -std::get<4>(k) } );
          if( it != chan.end() ) partner[c] = it->second;
        }
        auto cscale = [&](int c) { const int m = std::get<4>( chan_keys[c] ); return m < 0 && ( -m ) % 2 == 1 ? -1.0 : 1.0; };
        std::vector<bool> done( mono_keys.size(), false );
        struct Mono { std::vector<int> f; double c; };
        std::vector<Mono> monos;
        for(size_t k=0;k<mono_keys.size();k++)
        {
          if( done[k] ) continue;
          done[k] = true;
          std::vector<std::pair<int,int>> pk;
          double ratio = 1.0;
          bool ok = true;
          for(const auto& fp : mono_keys[k])
          {
            if( partner[fp.first] < 0 ) { ok = false; break; }
            pk.emplace_back( partner[fp.first], fp.second );
            for(int e=0;e<fp.second;e++) ratio *= cscale( fp.first ) / cscale( partner[fp.first] );
          }
          double c = mono_coef[k];
          if( ok )
          {
            std::sort( pk.begin(), pk.end() );
            const auto it = mono.find( pk );
            if( it != mono.end() && size_t(it->second) != k )
            {
              const double cp = mono_coef[ it->second ];
              const double sc = std::max( { 1.0, std::fabs(c), std::fabs(cp) } );
              if( std::fabs( c - ratio * cp ) > 1.0e-10 * sc ) fail(wb, ": the coefficients are not real-symmetric (complex density)");
              c += ratio * cp;
              done[ it->second ] = true;
            }
            else if( it == mono.end() && std::fabs( c ) > 0.0 ) fail(wb, ": a monomial without its complex-conjugate partner (complex density)");
          }
          Mono x;
          for(const auto& fp : mono_keys[k]) for(int e=0;e<fp.second;e++) x.f.push_back( fp.first );
          std::sort( x.f.begin(), x.f.end() );
          if( int(x.f.size()) > kYe3tMaxRank ) fail(wb, ": rank ", x.f.size(), " > ", kYe3tMaxRank);
          x.c = c;
          monos.push_back( std::move(x) );
        }
        // prefix tree of the factor lists (powers expanded), pre-order
        std::stable_sort( monos.begin(), monos.end(), [](const Mono& a, const Mono& b) { return a.f < b.f; } );
        std::vector<int> path;
        int nnodes = 0;
        for(const auto& x : monos)
        {
          size_t cp = 0;
          while( cp < path.size() && cp < x.f.size() && path[cp] == x.f[cp] ) ++ cp;
          if( cp == x.f.size() && cp == path.size() ) { node_coef.back() += x.c; continue; }
          path.resize( cp );
          for(size_t d=cp;d<x.f.size();d++)
          {
            path.push_back( x.f[d] );
            node_ch.push_back( x.f[d] );
            node_depth.push_back( d + 1 );
            node_coef.push_back( 0.0 );
            ++ nnodes;
          }
          node_coef.back() = x.c;
        }
        nmono += monos.size();
        nnode += nnodes;
        src0.push_back( src_kind.size() );
        ch0.push_back( ch_src.size() );
        no0.push_back( no0.back() + nnodes );
        T->nsrc_max = std::max( T->nsrc_max , int( skeys.size() ) );
        T->nch_max = std::max( T->nch_max , int( chan_keys.size() ) );
      }
      T->src_kind.assign( src_kind.begin(), src_kind.end() );
      T->src_r.assign( src_r.begin(), src_r.end() );
      T->src_lm.assign( src_lm.begin(), src_lm.end() );
      T->src_nu.assign( src_nu.begin(), src_nu.end() );
      T->ch_src.assign( ch_src.begin(), ch_src.end() );
      T->ch_tr.assign( ch_tr.begin(), ch_tr.end() );
      T->node_ch.assign( node_ch.begin(), node_ch.end() );
      T->node_depth.assign( node_depth.begin(), node_depth.end() );
      T->node_coef.assign( node_coef.begin(), node_coef.end() );
      T->src0.assign( src0.begin(), src0.end() );
      T->srcb.assign( srcb.begin(), srcb.end() );
      T->ch0.assign( ch0.begin(), ch0.end() );
      T->no0.assign( no0.begin(), no0.end() );

      lout << "YE3T linear ACE: " << nel << " element(s), " << nfunc << " functions (" << nrows << " rows -> " << nmono
           << " real monomials, " << nnode << " product nodes), lmax " << T->lmax << ", rcut " << T->rcut << " ang, " << src_kind.size() << " sources, "
           << T->bytes() / 1024 << " KB of tables" << ( T->analytic ? ", analytic radials" : "" ) << std::endl;
      return T;
    }

    // inverse of a small complex matrix (Gauss-Jordan, partial pivoting), row-major n x n
    std::vector< std::complex<double> > invert( std::vector< std::complex<double> > m, int n )
    {
      std::vector< std::complex<double> > inv( size_t(n) * n, 0.0 );
      for(int k=0;k<n;k++) inv[size_t(k)*n+k] = 1.0;
      for(int c=0;c<n;c++)
      {
        int piv = c;
        for(int r=c+1;r<n;r++) if( std::abs( m[size_t(r)*n+c] ) > std::abs( m[size_t(piv)*n+c] ) ) piv = r;
        if( std::abs( m[size_t(piv)*n+c] ) < 1e-14 ) fail("singular real-form matrix");
        for(int k=0;k<n;k++) { std::swap( m[size_t(c)*n+k], m[size_t(piv)*n+k] ); std::swap( inv[size_t(c)*n+k], inv[size_t(piv)*n+k] ); }
        const std::complex<double> d = 1.0 / m[size_t(c)*n+c];
        for(int k=0;k<n;k++) { m[size_t(c)*n+k] *= d; inv[size_t(c)*n+k] *= d; }
        for(int r=0;r<n;r++)
        {
          if( r == c ) continue;
          const std::complex<double> f = m[size_t(r)*n+c];
          if( f == 0.0 ) continue;
          for(int k=0;k<n;k++) { m[size_t(r)*n+k] -= f * m[size_t(c)*n+k]; inv[size_t(r)*n+k] -= f * inv[size_t(c)*n+k]; }
        }
      }
      return inv;
    }

    // tagged-Cauchy slice v2 (ye3t-lammps model_family tagged_cauchy): flat tables of ye3t_tagged.h. B: the .yace
    // backbone of a composite (same elements, cutoffs of the edge count), or null
    std::shared_ptr<Ye3tTaggedTables> load_tagged( const YAML::Node& root, const Ye3tTables* B )
    {
      const std::string w = "tagged component";
      if( req(root, "schema", w).as<std::string>() != "ye3t_tagged_cauchy_slice_v2" ) fail("tagged component: schema ye3t_tagged_cauchy_slice_v2 expected");
      auto TT = std::make_shared<Ye3tTaggedTables>();
      Ye3tTaggedTables& t = *TT;
      t.elements = req(root, "species_order", w).as< std::vector<std::string> >();
      const int nel = t.nel = t.elements.size();
      if( nel == 0 ) fail("tagged component: no species");
      if( B && B->elements != t.elements ) fail("ordinary and tagged components: different elements or element order");
      t.cutoff = req(root, "cutoff", w).as<double>();
      if( B && std::abs( B->rcut - t.cutoff ) > 2e-12 * std::max( 1.0 , t.cutoff ) ) fail("ordinary and tagged components: different cutoffs (", B->rcut, ", ", t.cutoff, ")");
      t.tag_count = req(root, "tag_count", w).as<int>();
      for(const char* k : { "combination_matrix", "descriptor_selection" }) if( root[k] && !root[k].IsNull() ) fail("tagged component: '", k, "' not supported");
      if( root["pooling"] && root["pooling"].as<std::string>() != "ordered_sum" ) fail("tagged component: pooling ", root["pooling"].as<std::string>(), " not supported");
      if( root["context_policy"] && root["context_policy"].as<std::string>() != "inclusive" ) fail("tagged component: context_policy ", root["context_policy"].as<std::string>(), " not supported");

      // radial basis
      const YAML::Node rd = req(root, "radial_definition", w);
      if( req(rd, "kind", "radial_definition").as<std::string>() != "pace_cheb_exp_cos" ) fail("tagged component: radial kind pace_cheb_exp_cos expected");
      const YAML::Node rp = req(rd, "parameters", "radial_definition");
      t.rc_rad = req(rp, "rc", "radial parameters").as<double>();
      t.rad.rc = t.rc_rad;
      t.rad.dcut = req(rp, "cutoff_width", "radial parameters").as<double>();
      t.rad.lambda = req(rp, "lmbda", "radial parameters").as<double>();
      t.rad.nbase = req(rp, "radial_count", "radial parameters").as<int>();
      if( t.rad.nbase < 1 ) fail("tagged component: radial_count < 1");

      // real forms: real coordinate a = sum_row U^-1[a][row] Y_{m(row)} (orthonormal Y) -> coefficients of Re / Im Y_l|m|
      // in the Y_00 = 1 normalisation (Y_{l,-m} = (-1)^m conj Y_lm)
      struct Form { int l = 0; std::vector<double> cf; };   // [2l+1][l+1][2]
      std::map<std::string, Form> forms;
      const double kappa = 1.0 / std::sqrt( 4.0 * M_PI );
      // real_form_records {id: record}, or (older slices) compiled_artifact.payload.real_forms [record with real_form_id]
      std::vector< std::pair<std::string, YAML::Node> > records;
      if( root["real_form_records"] ) for(const auto& kv : root["real_form_records"]) records.emplace_back( kv.first.as<std::string>(), kv.second );
      else if( root["compiled_artifact"] && root["compiled_artifact"]["payload"] && root["compiled_artifact"]["payload"]["real_forms"] )
        for(const auto& r : root["compiled_artifact"]["payload"]["real_forms"]) records.emplace_back( req(r, "real_form_id", "real form").as<std::string>(), r );
      else fail("tagged component: no real form (real_form_records or compiled_artifact.payload.real_forms)");
      for(const auto& kv : records)
      {
        const YAML::Node r = kv.second;
        Form f;
        f.l = req(r, "angular_l", "real form").as<int>();
        const int n = 2 * f.l + 1;
        const auto mo = req(r, "magnetic_order", "real form").as< std::vector<int> >();
        const YAML::Node U = req(r, "real_to_complex_matrix", "real form");
        if( int(mo.size()) != n || int(U.size()) != n ) fail("real form ", kv.first, ": bad size");
        std::vector< std::complex<double> > u( size_t(n) * n );
        for(int row=0;row<n;row++)
        {
          if( int(U[row].size()) != n ) fail("real form: bad row");
          for(int col=0;col<n;col++)
          {
            const auto b = req(U[row][col], "binary64", "real form").as< std::vector<double> >();
            u[size_t(row)*n+col] = std::complex<double>( b.at(0), b.at(1) );
          }
        }
        const auto ui = invert( u, n );
        f.cf.assign( size_t(n) * ( f.l + 1 ) * 2, 0.0 );
        for(int a=0;a<n;a++)
          for(int row=0;row<n;row++)
          {
            const std::complex<double> c = kappa * ui[size_t(a)*n+row];
            const int m = mo[row], am = std::abs( m );
            if( am > f.l ) fail("real form: bad magnetic order");
            double* o = f.cf.data() + ( size_t(a) * ( f.l + 1 ) + am ) * 2;
            if( m >= 0 ) { o[0] += c.real(); o[1] -= c.imag(); }
            else { const double sg = am % 2 ? -1.0 : 1.0; o[0] += sg * c.real(); o[1] += sg * c.imag(); }
          }
        forms[ kv.first ] = f;
      }

      // channels, components
      const YAML::Node chs = req(root, "channel_real_forms", w);
      const int nch = chs.size();
      std::vector<int> ch_off( nch + 1, -1 );
      std::vector<YAML::Node> by_index( nch );
      std::vector<bool> seen( nch, false );
      for(size_t q=0;q<chs.size();q++)
      {
        const YAML::Node c = chs[q];
        const int ci = req(c, "channel_index", "channel").as<int>();
        if( ci < 0 || ci >= nch || seen[ci] ) fail("channels: bad or duplicate channel_index");
        seen[ci] = true;
        by_index[ci].reset( c );
      }
      std::vector<int> comp_rad, comp_nu, comp_l, comp_cf;
      std::vector<double> cf;
      for(int ci=0;ci<nch;ci++)
      {
        const YAML::Node c = by_index[ci];
        const int l = req(c, "l", "channel").as<int>(), n = req(c, "radial_channel", "channel").as<int>();
        const std::string nus = req(c, "neighbor_species", "channel").as<std::string>();
        const auto it = std::find( t.elements.begin(), t.elements.end(), nus );
        if( it == t.elements.end() ) fail("channel: unknown neighbour species ", nus);
        const auto fi = forms.find( req(c, "real_form_id", "channel").as<std::string>() );
        if( fi == forms.end() || fi->second.l != l ) fail("channel ", ci, ": real form missing or of another l");
        if( n < 0 || n >= t.rad.nbase ) fail("channel ", ci, ": radial channel out of range");
        if( l > kYe3tMaxLmax ) fail("channel ", ci, ": l > ", kYe3tMaxLmax);
        t.lmax = std::max( t.lmax , l );
        ch_off[ci] = comp_rad.size();
        for(int a=0;a<2*l+1;a++)
        {
          comp_rad.push_back( n ); comp_nu.push_back( it - t.elements.begin() ); comp_l.push_back( l ); comp_cf.push_back( cf.size() );
          cf.insert( cf.end(), fi->second.cf.begin() + size_t(a) * ( l + 1 ) * 2, fi->second.cf.begin() + size_t(a+1) * ( l + 1 ) * 2 );
        }
      }
      ch_off[nch] = comp_rad.size();
      t.ncomp = comp_rad.size();
      auto flat = [&](const YAML::Node& key) -> int
      {
        const auto k = key.as< std::vector<int> >();
        if( k.size() != 2 || k[0] < 0 || k[0] >= nch || k[1] < 0 || k[1] >= ch_off[k[0]+1] - ch_off[k[0]] ) fail("bad (channel, component) key");
        return ch_off[k[0]] + k[1];
      };

      // densities, moments (prefix tree)
      const YAML::Node prog = req(root, "real_moment_program", w);
      if( prog["tag_count"] && prog["tag_count"].as<int>() != t.tag_count ) fail("real_moment_program: tag_count differs");
      std::vector<int> den_comp;
      for(const auto& k : req(prog, "real_density_keys", "real_moment_program")) den_comp.push_back( flat( k ) );
      t.nden = den_comp.size();
      std::map< std::vector<int>, std::vector<int> > mkeys;
      std::vector<int> mzero;
      {
        int m = 0;
        for(const auto& key : req(prog, "real_moment_keys", "real_moment_program"))
        {
          std::vector<int> f;
          for(const auto& k : key) f.push_back( flat( k ) );
          std::sort( f.begin(), f.end() );
          if( int(f.size()) > kTaggedMaxDepth ) fail("moment key of ", f.size(), " factors (", kTaggedMaxDepth, " supported)");
          if( f.empty() ) mzero.push_back( m ); else mkeys[f].push_back( m );
          ++ m;
        }
        t.nmom = m;
      }
      t.nvar = t.nden + t.nmom;
      std::vector<int> mnode_comp, mnode_depth, mt0, mterm;
      {
        std::vector<int> path;
        for(const auto& kv : mkeys)
        {
          const auto& f = kv.first;
          size_t c = 0;
          while( c < path.size() && c < f.size() && path[c] == f[c] ) ++ c;
          if( c == f.size() ) fail("moment tree: duplicate key");   // map keys are distinct, a prefix comes first
          for(size_t d=c;d<f.size();d++) { mnode_comp.push_back( f[d] ); mnode_depth.push_back( d + 1 ); mt0.push_back( mterm.size() ); }
          mterm.insert( mterm.end(), kv.second.begin(), kv.second.end() );
          path = f;
        }
        mt0.push_back( mterm.size() );
      }

      // readout: beta per element, terms -> (p, beta coef) per sorted variable list
      const YAML::Node bt = req(root, "beta", w);
      std::vector< std::vector<double> > beta( nel );
      if( bt.IsSequence() ) { if( nel != 1 ) fail("beta: a list is only accepted for one species"); beta[0] = bt.as< std::vector<double> >(); }
      else for(int e=0;e<nel;e++) beta[e] = req(bt, t.elements[e].c_str(), "beta").as< std::vector<double> >();
      std::vector<double> offset( nel, 0.0 );
      if( root["offsets"] && !root["offsets"].IsNull() ) for(int e=0;e<nel;e++) if( root["offsets"][t.elements[e]] ) offset[e] = root["offsets"][t.elements[e]].as<double>();
      std::vector< std::map< std::vector<int>, std::map<int,double> > > rterms( nel );
      size_t nterm = 0;
      for(const auto& term : req(prog, "terms", "real_moment_program"))
      {
        const int fe = req(term, "feature_index", "term").as<int>(), p = req(term, "p", "term").as<int>();
        const double coef = req(term, "coefficient", "term").as<double>();
        if( p < 0 || p > t.tag_count ) fail("term: p out of range");
        std::vector<int> v;
        for(int k : req(term, "density_factor_indices", "term").as< std::vector<int> >()) { if( k < 0 || k >= t.nden ) fail("term: density index out of range"); v.push_back( k ); }
        for(int k : req(term, "moment_indices", "term").as< std::vector<int> >()) { if( k < 0 || k >= t.nmom ) fail("term: moment index out of range"); v.push_back( t.nden + k ); }
        std::sort( v.begin(), v.end() );
        if( int(v.size()) > kTaggedMaxDepth ) fail("term of ", v.size(), " factors (", kTaggedMaxDepth, " supported)");
        for(int e=0;e<nel;e++)
        {
          if( fe < 0 || fe >= int( beta[e].size() ) ) fail("term: feature_index out of range of beta");
          rterms[e][v][p] += beta[e][fe] * coef;
        }
        t.pmax = std::max( t.pmax , p );
        ++ nterm;
      }
      std::vector<int> ro0{ 0 }, rnode_var, rnode_depth, rl0, rl_p, r0terms{ 0 }, r0t_p;
      std::vector<double> rl_coef, r0t_coef;
      for(int e=0;e<nel;e++)
      {
        std::vector<int> path;
        for(const auto& kv : rterms[e])
        {
          const auto& v = kv.first;
          if( v.empty() ) { for(const auto& pc : kv.second) { r0t_p.push_back( pc.first ); r0t_coef.push_back( pc.second ); } continue; }
          size_t c = 0;
          while( c < path.size() && c < v.size() && path[c] == v[c] ) ++ c;
          for(size_t d=c;d<v.size();d++) { rnode_var.push_back( v[d] ); rnode_depth.push_back( d + 1 ); rl0.push_back( rl_p.size() ); }
          for(const auto& pc : kv.second) { rl_p.push_back( pc.first ); rl_coef.push_back( pc.second ); }
          path = v;
        }
        ro0.push_back( rnode_var.size() );
        r0terms.push_back( r0t_p.size() );
      }
      rl0.push_back( rl_p.size() );

      // edge count cutoffs
      std::vector<double> zc2( size_t(nel) * nel );
      for(int a=0;a<nel;a++) for(int b=0;b<nel;b++)
      {
        const double c = B ? std::max( t.cutoff , B->bonds[size_t(a)*nel+b].rc ) : t.cutoff;
        zc2[size_t(a)*nel+b] = c * c;
      }

      t.nlh = ( t.lmax + 1 ) * ( t.lmax + 2 ) / 2;
      t.ylm_a.resize( t.nlh ); t.ylm_b.resize( t.nlh );
      ylm_coefficients( t.lmax, t.ylm_a.data(), t.ylm_b.data() );
      t.zcut2.assign( zc2.begin(), zc2.end() ); t.offset.assign( offset.begin(), offset.end() );
      t.comp_rad.assign( comp_rad.begin(), comp_rad.end() ); t.comp_nu.assign( comp_nu.begin(), comp_nu.end() );
      t.comp_l.assign( comp_l.begin(), comp_l.end() ); t.comp_cf.assign( comp_cf.begin(), comp_cf.end() ); t.cf.assign( cf.begin(), cf.end() );
      t.den_comp.assign( den_comp.begin(), den_comp.end() );
      t.mnode_comp.assign( mnode_comp.begin(), mnode_comp.end() ); t.mnode_depth.assign( mnode_depth.begin(), mnode_depth.end() );
      t.mt0.assign( mt0.begin(), mt0.end() ); t.mterm.assign( mterm.begin(), mterm.end() ); t.mzero.assign( mzero.begin(), mzero.end() );
      {
        std::vector<int> par( mnode_comp.size(), -1 ), last( kTaggedMaxDepth + 1, -1 );
        for(size_t n=0;n<mnode_comp.size();n++) { const int d = mnode_depth[n]; par[n] = d > 1 ? last[d-1] : -1; last[d] = n; }
        t.mnode_parent.assign( par.begin(), par.end() );
        // reverse sweep: the first child met (largest index) of each parent initialises it
        std::vector<int> first( par.size(), 0 ), hasch( par.size(), 0 ), seen( par.size(), 0 );
        for(int n=int(par.size())-1;n>=0;n--) if( par[n] >= 0 ) { hasch[par[n]] = 1; if( !seen[par[n]] ) { seen[par[n]] = 1; first[n] = 1; } }
        t.mnode_first.assign( first.begin(), first.end() );
        t.mnode_haschild.assign( hasch.begin(), hasch.end() );
        std::vector<int> slot( par.size(), -1 ), pslot( par.size(), -1 );
        int ns = 0;
        for(size_t n=0;n<par.size();n++) if( hasch[n] ) slot[n] = ns++;
        for(size_t n=0;n<par.size();n++) pslot[n] = par[n] >= 0 ? slot[par[n]] : -1;
        t.nslot = ns;
        t.mnode_slot.assign( slot.begin(), slot.end() );
        t.mnode_pslot.assign( pslot.begin(), pslot.end() );
      }
      {
        // distinct angular forms of the components
        std::map< std::pair<int, std::vector<double> >, int > amap;
        std::vector<int> al, acf, ca;
        for(int c=0;c<t.ncomp;c++)
        {
          const int l = comp_l[c];
          std::vector<double> row( cf.begin() + comp_cf[c], cf.begin() + comp_cf[c] + 2 * ( l + 1 ) );
          auto it = amap.find( { l, row } );
          if( it == amap.end() ) { it = amap.emplace( std::make_pair( l, row ), int( al.size() ) ).first; al.push_back( l ); acf.push_back( comp_cf[c] ); }
          ca.push_back( it->second );
        }
        t.nang = al.size();
        t.ang_l.assign( al.begin(), al.end() ); t.ang_cf.assign( acf.begin(), acf.end() ); t.comp_ang.assign( ca.begin(), ca.end() );
      }
      t.ro0.assign( ro0.begin(), ro0.end() ); t.rnode_var.assign( rnode_var.begin(), rnode_var.end() ); t.rnode_depth.assign( rnode_depth.begin(), rnode_depth.end() );
      t.rl0.assign( rl0.begin(), rl0.end() ); t.rl_p.assign( rl_p.begin(), rl_p.end() ); t.rl_coef.assign( rl_coef.begin(), rl_coef.end() );
      t.r0terms.assign( r0terms.begin(), r0terms.end() ); t.r0t_p.assign( r0t_p.begin(), r0t_p.end() ); t.r0t_coef.assign( r0t_coef.begin(), r0t_coef.end() );

      lout << "YE3T tagged-Cauchy (slice v2): " << nel << " element(s), tag_count " << t.tag_count << ", " << nch << " channels (" << t.ncomp
           << " components, lmax " << t.lmax << ", " << t.rad.nbase << " radials, rc " << t.rc_rad << " ang), " << t.nden << " densities, " << t.nmom
           << " moments (" << mnode_comp.size() << " product nodes), " << nterm << " terms (" << rnode_var.size() << " readout nodes), "
           << TT->bytes() / 1024 << " KB of tables" << std::endl;
      return TT;
    }

    inline void execute() override final
    {
      const std::string path = onika::data_file_path( parameters->file );
      std::shared_ptr<Ye3tTables> T;
      std::shared_ptr<Ye3tTaggedTables> TT;
      std::vector<std::string> elements;
      double rcut = 0.0;
      if( path.size() > 5 && path.substr( path.size() - 5 ) == ".json" )
      {
        m_fname = path;
        YAML::Node root;
        try { root = YAML::LoadFile( path ); }
        catch( const YAML::Exception& e ) { fail("cannot read model file: ", e.what()); }
        const std::string schema = req(root, "schema", "model").as<std::string>();
        const std::string dir = path.find('/') != std::string::npos ? path.substr( 0 , path.rfind('/') + 1 ) : std::string();
        std::string tagged_path = path;
        YAML::Node troot = root;
        if( schema == "ye3t_tagged_cauchy_composite_v1" )
        {
          const std::string bp = dir + req(req(root, "ordinary_component", "model"), "path", "ordinary_component").as<std::string>();
          tagged_path = dir + req(req(root, "tagged_component", "model"), "path", "tagged_component").as<std::string>();
          T = load_yace( bp );
          m_fname = tagged_path;
          try { troot = YAML::LoadFile( tagged_path ); }
          catch( const YAML::Exception& e ) { fail("cannot read tagged component: ", e.what()); }
        }
        else if( schema != "ye3t_tagged_cauchy_slice_v2" ) fail("unsupported schema ", schema, " (ye3t_tagged_cauchy_composite_v1 or ye3t_tagged_cauchy_slice_v2)");
        TT = load_tagged( troot, T.get() );
        elements = TT->elements;
        rcut = TT->cutoff;
        if( T ) rcut = std::max( rcut , T->rcut );
      }
      else
      {
        T = load_yace( path );
        elements = T->elements;
        rcut = T->rcut;
      }

      // exaStamp species -> model elements, by name
      Ye3tContext& ctx = *ye3t_ctx;
      ctx = Ye3tContext{};
      const auto& sp = *species;
      ctx.species_type.assign( sp.size(), -1 );
      for(size_t s=0;s<sp.size();s++)
      {
        const auto it = std::find( elements.begin(), elements.end(), std::string( sp[s].m_name ) );
        if( it != elements.end() ) ctx.species_type[s] = it - elements.begin();
        else lout << "ye3t_init: warning: species " << sp[s].m_name << " is not an element of the model, its atoms are ignored" << std::endl;
      }
      ctx.tables = T;
      ctx.tagged = TT;
      ctx.rcut = rcut;
      ctx.mem_gb = parameters->mem_gb;
      ctx.emulate = parameters->emulate;
      *rcut_max = std::max( *rcut_max , rcut );
    }
  };

  ONIKA_AUTORUN_INIT(ye3t_init)
  {
    OperatorNodeFactory::instance()->register_factory( "ye3t_init", make_simple_operator<Ye3tInit> );
  }

}
