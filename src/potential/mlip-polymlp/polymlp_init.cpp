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
#include <exaStamp/mlip_polymlp/polymlp_context.h>
#include <exaStamp/mlip_utils/mlip_ylm.h>

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
  using namespace mlip_polymlp;

  class PolymlpInit : public OperatorNode
  {
    ADD_SLOT( PolymlpParams   , parameters  , INPUT        , REQUIRED , DocString{"file: model exported by scripts/polymlp_export; mem_gb, emulate: see documentation"} );
    ADD_SLOT( double          , rcut_max    , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( ParticleSpecies , species     , INPUT        , REQUIRED );
    ADD_SLOT( PolymlpContext  , polymlp_ctx , OUTPUT );

    std::string m_fname;

    template<class... T> [[noreturn]] void fail(const T&... msg)
    {
      std::ostringstream oss;
      ( oss << ... << msg );
      fatal_error() << "polymlp_init: " << oss.str() << " (" << m_fname << ")" << std::endl;
      std::abort();
    }

    YAML::Node req(const YAML::Node& n, const char* key)
    {
      const YAML::Node x = n[key];
      if( !x ) fail("missing '", key, "'");
      return x;
    }

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(

Polynomial machine-learning potential (A. Seko; models of pypolymlp and of the polynomial MLP repository, used by
lammps-polymlp-package as pair_style polymlp), for polymlp_force. Native onika kernels (GPU if available, else OpenMP).
The polymlp file (polymlp.lammps or polymlp.yaml) is first converted by the exaStamp tool scripts/polymlp_export,
which enumerates the invariant features with the polymlp library:

  polymlp_export polymlp.lammps polymlp.exastamp.yaml

polymlp_init:
  parameters:
    file: polymlp.exastamp.yaml
    mem_gb: 1.0       # budget of the per-chunk buffers
    emulate: 0        # debug: N > 0 runs the GPU code paths on the host with blocks of N threads (slow)

The exaStamp species are mapped to the model elements by name. Gaussian radial functions only.

)EOF";
    }

    inline void execute() override final
    {
      m_fname = onika::data_file_path( parameters->file );
      YAML::Node root;
      try { root = YAML::LoadFile( m_fname ); }
      catch( const YAML::Exception& e ) { fail("cannot read model file: ", e.what()); }
      if( !root["format"] || root["format"].as<std::string>() != "exastamp-polymlp" ) fail("not a model written by scripts/polymlp_export");
      if( root["version"].as<int>() != 1 ) fail("unsupported version ", root["version"].as<int>());

      auto T = std::make_shared<PolymlpTables>();
      const auto elements = req(root, "elements").as< std::vector<std::string> >();
      const int nel = elements.size();
      if( nel == 0 ) fail("no element");
      T->elements = elements;
      T->nel = nel;
      T->rc = req(root, "cutoff").as<double>();
      const std::string ftype = req(root, "feature_type").as<std::string>();
      if( ftype != "gtinv" && ftype != "pair" ) fail("unsupported feature_type ", ftype);
      T->gtinv = ftype == "gtinv";
      const auto gs = req(root, "gaussians").as< std::vector< std::vector<double> > >();
      T->ng = gs.size();
      if( T->ng == 0 || T->ng > kPolymlpMaxRadial ) fail(T->ng, " gaussians, 1 to ", kPolymlpMaxRadial, " supported");
      for(const auto& g : gs) { if( g.size() != 2 ) fail("gaussians: [beta, mu] expected"); T->beta.push_back( g[0] ); T->mu.push_back( g[1] ); }
      std::vector< std::pair<int,int> > tp_types;
      for(const auto& tp : req(root, "type_pairs"))
      {
        const auto ty = tp["types"].as< std::vector<int> >();
        if( ty.size() != 2 || ty[0] < 0 || ty[0] >= nel || ty[1] < 0 || ty[1] >= nel ) fail("type_pairs: bad types");
        tp_types.emplace_back( ty[0], ty[1] );
      }

      const YAML::Node types = req(root, "types");
      if( int(types.size()) != nel ) fail("expected ", nel, " types");
      const double kappa = 1.0 / std::sqrt( 4.0 * M_PI );
      std::vector<int> src0{0}, srcb, ch0{0}, ft0{0}, po0{0}, no0{0}, src_n, src_lm, src_nu, ch_src, ch_conj, node_ch, node_depth, nt_feat, pfac;
      std::vector<double> ch_scale, nt_coef, pcoef;
      std::vector<int64_t> nto{0}, pfo{0};
      size_t nfeat = 0, nterm = 0, npoly = 0, nnode = 0;
      for(int t=0;t<nel;t++)
      {
        const YAML::Node ty = types[t];
        const int w = T->gtinv ? 5 : 2;
        const auto op = req(ty, "order_parameters").as< std::vector<int> >();
        if( op.size() % w != 0 ) fail("type ", t, ": order_parameters size");
        const int nop = op.size() / w;
        // sources (nu, n, l, |m|), sorted by neighbour type, in order of first use
        using SKey = std::tuple<int,int,int,int>;
        std::vector<SKey> skeys;
        std::map<SKey,int> smap;
        std::vector<SKey> okeys( nop );
        std::vector<double> oscale( nop );
        std::vector<int> oconj( nop );
        for(int k=0;k<nop;k++)
        {
          const int* r = op.data() + size_t(k) * w;
          const int n = r[0], tp = T->gtinv ? r[3] : r[1];
          const int l = T->gtinv ? r[1] : 0, m = T->gtinv ? r[2] : 0;
          if( n < 0 || n >= T->ng ) fail("type ", t, ": gaussian index out of range");
          if( tp < 0 || tp >= int(tp_types.size()) ) fail("type ", t, ": type pair out of range");
          const auto [ta, tb] = tp_types[tp];
          if( ta != t && tb != t ) fail("type ", t, ": order parameter of a type pair without this type");
          const int nu = ta == t ? tb : ta;
          if( l < 0 || l > kPolymlpMaxLmax || std::abs(m) > l ) fail("type ", t, ": l = ", l, " m = ", m, " (lmax ", kPolymlpMaxLmax, ")");
          T->lmax = std::max( T->lmax , l );
          okeys[k] = SKey{ nu, n, l, std::abs(m) };
          if( smap.find( okeys[k] ) == smap.end() ) { smap.emplace( okeys[k], 0 ); skeys.push_back( okeys[k] ); }
          // polymlp: orthonormal harmonics of r_i - r_j: kappa (-1)^l Y_lm(r_ij); negative m: (-1)^m conj
          double s = T->gtinv ? kappa * ( l % 2 ? -1.0 : 1.0 ) : 1.0;
          if( m < 0 && ( -m ) % 2 == 1 ) s = -s;
          oscale[k] = s;
          oconj[k] = m < 0 ? 1 : 0;
        }
        std::stable_sort( skeys.begin(), skeys.end(), [](const SKey& a, const SKey& b) { return std::get<0>(a) < std::get<0>(b); } );
        for(size_t q=0;q<skeys.size();q++) smap[ skeys[q] ] = q;
        for(int nu=0;nu<=nel;nu++)
        {
          int q = 0;
          while( q < int(skeys.size()) && std::get<0>( skeys[q] ) < nu ) ++ q;
          srcb.push_back( q );
        }
        for(const auto& s : skeys) { src_nu.push_back( std::get<0>(s) ); src_n.push_back( std::get<1>(s) ); src_lm.push_back( mlip_utils::lmh( std::get<2>(s), std::get<3>(s) ) ); }
        for(int k=0;k<nop;k++) { ch_src.push_back( smap[ okeys[k] ] ); ch_conj.push_back( oconj[k] ); ch_scale.push_back( oscale[k] ); }

        // features: terms sorted by factor list -> prefix tree in pre-order (shared products)
        const auto fo = req(ty, "feature_offsets").as< std::vector<int64_t> >();
        const auto tc = req(ty, "term_coefficients").as< std::vector<double> >();
        const auto to = req(ty, "term_offsets").as< std::vector<int64_t> >();
        const auto ti = req(ty, "term_ids").as< std::vector<int> >();
        if( fo.empty() || fo.front() != 0 || size_t(fo.back()) != tc.size() || to.size() != tc.size() + 1 || size_t(to.back()) != ti.size() ) fail("type ", t, ": inconsistent feature tables");
        const int nf = fo.size() - 1;
        // complex-conjugate partner of each channel (same source, other conjugation; m = 0: itself): the term with the
        // partner channels has the conjugate product times prod(s_c / s_partner), the same real part -> one canonical term
        std::map< std::pair<SKey,int>, int > chan_of;
        for(int k=0;k<nop;k++) chan_of[ { okeys[k], oconj[k] } ] = k;
        std::vector<int> partner( nop, -1 );
        for(int k=0;k<nop;k++)
        {
          if( std::get<3>( okeys[k] ) == 0 ) { partner[k] = k; continue; }
          const auto it = chan_of.find( { okeys[k], 1 - oconj[k] } );
          if( it != chan_of.end() ) partner[k] = it->second;
        }
        struct Term { std::vector<int> f; int k; double c; };
        std::vector<Term> terms;
        terms.reserve( tc.size() );
        for(int k=0;k<nf;k++) for(int64_t s=fo[k];s<fo[k+1];s++)
        {
          Term x{ std::vector<int>( ti.begin() + to[s], ti.begin() + to[s+1] ), k, tc[s] };
          if( x.f.empty() || int(x.f.size()) > kPolymlpMaxFactors ) fail("type ", t, ": feature term with ", x.f.size(), " factors (1 to ", kPolymlpMaxFactors, ")");
          for(int id : x.f) if( id < 0 || id >= nop ) fail("type ", t, ": feature term id out of range");
          std::sort( x.f.begin(), x.f.end() );
          std::vector<int> g;
          double ratio = 1.0;
          bool ok = true;
          for(int id : x.f) { if( partner[id] < 0 ) { ok = false; break; } g.push_back( partner[id] ); ratio *= oscale[id] / oscale[partner[id]]; }
          if( ok )
          {
            std::sort( g.begin(), g.end() );
            if( g < x.f ) { x.f = g; x.c *= ratio; }
          }
          terms.push_back( std::move(x) );
        }
        std::stable_sort( terms.begin(), terms.end(), [](const Term& x, const Term& y) { return x.f != y.f ? x.f < y.f : x.k < y.k; } );
        {
          // merge equal (factors, feature) terms
          size_t w = 0;
          for(size_t r=0;r<terms.size();r++)
          {
            if( w > 0 && terms[w-1].f == terms[r].f && terms[w-1].k == terms[r].k ) { terms[w-1].c += terms[r].c; continue; }
            if( w != r ) terms[w] = std::move( terms[r] );
            ++ w;
          }
          terms.resize( w );
        }
        std::vector<int> path;
        int nnodes = 0;
        for(const auto& x : terms)
        {
          size_t cp = 0;
          while( cp < path.size() && cp < x.f.size() && path[cp] == x.f[cp] ) ++ cp;
          if( cp == x.f.size() && cp == path.size() ) { nt_feat.push_back( x.k ); nt_coef.push_back( x.c ); nto.back() = nt_feat.size(); continue; }
          path.resize( cp );
          for(size_t d=cp;d<x.f.size();d++)
          {
            path.push_back( x.f[d] );
            node_ch.push_back( x.f[d] );
            node_depth.push_back( d + 1 );
            nto.push_back( nt_feat.size() );
            ++ nnodes;
          }
          nt_feat.push_back( x.k ); nt_coef.push_back( x.c ); nto.back() = nt_feat.size();
        }
        no0.push_back( no0.back() + nnodes );
        nterm += terms.size(); nnode += nnodes;
        // polynomial
        const auto pc = req(ty, "poly_coefficients").as< std::vector<double> >();
        const auto po = req(ty, "poly_offsets").as< std::vector<int64_t> >();
        const auto pf = req(ty, "poly_features").as< std::vector<int> >();
        if( po.size() != pc.size() + 1 || size_t(po.back()) != pf.size() ) fail("type ", t, ": inconsistent polynomial tables");
        const int64_t pbase = pfac.size();
        for(size_t u=0;u<pc.size();u++)
        {
          if( po[u+1] - po[u] > kPolymlpMaxFactors ) fail("type ", t, ": polynomial term with more than ", kPolymlpMaxFactors, " factors");
          pcoef.push_back( pc[u] );
          pfo.push_back( pbase + po[u+1] );
        }
        for(int k : pf) { if( k < 0 || k >= nf ) fail("type ", t, ": polynomial feature out of range"); pfac.push_back( k ); }

        src0.push_back( src_n.size() );
        ch0.push_back( ch_src.size() );
        ft0.push_back( ft0.back() + nf );
        po0.push_back( pcoef.size() );
        T->nsrc_max = std::max( T->nsrc_max , int( skeys.size() ) );
        T->nfeat_max = std::max( T->nfeat_max , nf );
        T->nch_max = std::max( T->nch_max , nop );
        nfeat += nf; npoly += pc.size();
      }
      T->nlh = ( T->lmax + 1 ) * ( T->lmax + 2 ) / 2;
      T->ylm_a.resize( T->nlh );
      T->ylm_b.resize( T->nlh );
      mlip_utils::ylm_coefficients( T->lmax, T->ylm_a.data(), T->ylm_b.data() );
      T->src0.assign( src0.begin(), src0.end() ); T->srcb.assign( srcb.begin(), srcb.end() ); T->ch0.assign( ch0.begin(), ch0.end() );
      T->ft0.assign( ft0.begin(), ft0.end() ); T->po0.assign( po0.begin(), po0.end() );
      T->src_n.assign( src_n.begin(), src_n.end() ); T->src_lm.assign( src_lm.begin(), src_lm.end() ); T->src_nu.assign( src_nu.begin(), src_nu.end() );
      T->ch_src.assign( ch_src.begin(), ch_src.end() ); T->ch_conj.assign( ch_conj.begin(), ch_conj.end() ); T->ch_scale.assign( ch_scale.begin(), ch_scale.end() );
      T->no0.assign( no0.begin(), no0.end() ); T->node_ch.assign( node_ch.begin(), node_ch.end() ); T->node_depth.assign( node_depth.begin(), node_depth.end() );
      T->nto.assign( nto.begin(), nto.end() ); T->nt_feat.assign( nt_feat.begin(), nt_feat.end() ); T->nt_coef.assign( nt_coef.begin(), nt_coef.end() );
      T->pcoef.assign( pcoef.begin(), pcoef.end() ); T->pfo.assign( pfo.begin(), pfo.end() ); T->pfac.assign( pfac.begin(), pfac.end() );

      PolymlpContext& ctx = *polymlp_ctx;
      ctx = PolymlpContext{};
      const auto& sp = *species;
      ctx.species_type.assign( sp.size(), -1 );
      for(size_t s=0;s<sp.size();s++)
      {
        const auto it = std::find( elements.begin(), elements.end(), std::string( sp[s].m_name ) );
        if( it != elements.end() ) ctx.species_type[s] = it - elements.begin();
        else lout << "polymlp_init: warning: species " << sp[s].m_name << " is not an element of the model, its atoms are ignored" << std::endl;
      }
      ctx.tables = T;
      ctx.mem_gb = parameters->mem_gb;
      ctx.emulate = parameters->emulate;
      *rcut_max = std::max( *rcut_max , T->rc );
      lout << "polymlp: " << nel << " element(s), " << ftype << ", " << T->ng << " gaussians, lmax " << T->lmax << ", rcut " << T->rc << " ang, "
           << nfeat << " features (" << nterm << " terms, " << nnode << " product nodes), " << npoly << " polynomial terms, " << src_n.size() << " sources, "
           << T->bytes() / 1024 << " KB of tables" << std::endl;
    }
  };

  ONIKA_AUTORUN_INIT(polymlp_init)
  {
    OperatorNodeFactory::instance()->register_factory( "polymlp_init", make_simple_operator<PolymlpInit> );
  }

}
