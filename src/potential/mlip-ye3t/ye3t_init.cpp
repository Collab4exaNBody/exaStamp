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
    ADD_SLOT( Ye3tParams      , parameters , INPUT        , REQUIRED , DocString{"file: .yace model (linear ACE, ye3t-lammps pair_style ye3t); analytic, mem_gb, emulate: see documentation"} );
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

The exaStamp species are mapped to the model elements by name.

)EOF";
    }

    inline void execute() override final
    {
      m_fname = onika::data_file_path( parameters->file );
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
        std::vector<double> gv( size_t(b.nint) * b.nbase ), gd( gv.size() ), cv( size_t(b.nint) * b.ncontr ), cd( cv.size() );
        for(int k=0;k<b.nint;k++)
        {
          const double r = h * double( k + 1 );
          double* g = gv.data() + size_t(k) * b.nbase;
          double* dg = gd.data() + size_t(k) * b.nbase;
          if( r >= b.rc ) { std::fill( g, g + b.nbase, 0.0 ); std::fill( dg, dg + b.nbase, 0.0 ); }
          else cheb_exp_cos( b, r, g, dg );
          for(int q=0;q<b.ncontr;q++)
          {
            double s = 0.0, ds = 0.0;
            for(int kk=0;kk<b.nbase;kk++) { s += pool[ b.crad + size_t(q)*b.nbase + kk ] * g[kk]; ds += pool[ b.crad + size_t(q)*b.nbase + kk ] * dg[kk]; }
            cv[ size_t(k)*b.ncontr + q ] = s;
            cd[ size_t(k)*b.ncontr + q ] = ds;
          }
        }
        b.bspl = pool.size();
        pool.resize( pool.size() + size_t(b.nint + 1) * b.nbase * 4 );
        build_spline( gv, gd, b.nint, b.nbase, h, pool.data() + b.bspl );
        b.cspl = pool.size();
        pool.resize( pool.size() + size_t(b.nint + 1) * b.ncontr * 4 );
        build_spline( cv, cd, b.nint, b.ncontr, h, pool.data() + b.cspl );
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
      ctx.mem_gb = parameters->mem_gb;
      ctx.emulate = parameters->emulate;
      *rcut_max = std::max( *rcut_max , T->rcut );

      lout << "YE3T linear ACE: " << nel << " element(s), " << nfunc << " functions (" << nrows << " rows -> " << nmono
           << " real monomials, " << nnode << " product nodes), lmax " << T->lmax << ", rcut " << T->rcut << " ang, " << src_kind.size() << " sources, "
           << T->bytes() / 1024 << " KB of tables" << ( T->analytic ? ", analytic radials" : "" ) << std::endl;
    }
  };

  ONIKA_AUTORUN_INIT(ye3t_init)
  {
    OperatorNodeFactory::instance()->register_factory( "ye3t_init", make_simple_operator<Ye3tInit> );
  }

}
