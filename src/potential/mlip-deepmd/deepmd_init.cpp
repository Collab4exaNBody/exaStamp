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

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/file_utils.h>

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "deepmd_model.h"
#include "deepmd_dpa3_driver.h"
#include "deepmd_dpa2.h"
#include "deepmd_dpa4.h"

namespace exaStamp
{

  using namespace exanb;

  class DeepmdInit : public OperatorNode
  {
    ADD_SLOT( DeepmdParams    , parameters , INPUT        , REQUIRED , DocString{"exported DeePMD model file"} );
    ADD_SLOT( double          , rcut_max   , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( ParticleSpecies , species    , INPUT        , REQUIRED );
    ADD_SLOT( DeepmdContext   , deepmd_ctx , OUTPUT );

    std::string m_fname;

    template<class... T> [[noreturn]] void fail(const T&... msg)
    {
      std::ostringstream oss;
      ( oss << ... << msg );
      fatal_error() << "deepmd_init: " << oss.str() << " (" << m_fname << ")" << std::endl;
      std::abort();
    }

    static int activation_code(const std::string& a)
    {
      if( a == "tanh" ) return DP_ACT_TANH;
      if( a == "relu" ) return DP_ACT_RELU;
      if( a == "gelu" || a == "gelu_tf" ) return DP_ACT_GELU;
      if( a == "softplus" ) return DP_ACT_SOFTPLUS;
      if( a == "sigmoid" ) return DP_ACT_SIGMOID;
      if( a == "silu" ) return DP_ACT_SILU;
      if( a == "none" || a == "linear" ) return DP_ACT_NONE;
      return -1;
    }

    // appends the layers of one network to the context, returns its descriptor
    DPNet read_net(const YAML::Node& node, DeepmdContext& ctx, const std::string& what, int max_layers, int max_width, const char* limits)
    {
      const YAML::Node layers = node["layers"];
      if( !layers || !layers.IsSequence() || layers.size() == 0 ) fail(what, ": no layers");
      if( int(layers.size()) > max_layers ) fail(what, ": ", layers.size(), " layers, at most ", max_layers, " (", limits, "_LAYERS)");
      DPNet net;
      net.first_layer = ctx.layers.size();
      net.nlayers = layers.size();
      for(const auto& ln : layers)
      {
        DPLayer L;
        L.nin = ln["nin"].as<int>();
        L.nout = ln["nout"].as<int>();
        L.activation = activation_code( ln["activation"].as<std::string>() );
        if( L.activation < 0 ) fail(what, ": unsupported activation '", ln["activation"].as<std::string>(), "'");
        if( ln["resnet"].as<bool>() ) L.resnet = ( L.nout == L.nin ) ? 1 : ( L.nout == 2 * L.nin ) ? 2 : 0;
        if( L.nout > max_width ) fail(what, ": layer width ", L.nout, " exceeds ", max_width, " (", limits, "_WIDTH)");
        const auto w = ln["w"].as< std::vector<double> >();
        const auto b = ln["b"].as< std::vector<double> >();
        if( w.size() != size_t(L.nin) * L.nout || b.size() != size_t(L.nout) ) fail(what, ": inconsistent layer sizes");
        L.w = ctx.weights.size(); ctx.weights.insert( ctx.weights.end(), w.begin(), w.end() );
        L.b = ctx.weights.size(); ctx.weights.insert( ctx.weights.end(), b.begin(), b.end() );
        if( ln["idt"] )
        {
          const auto idt = ln["idt"].as< std::vector<double> >();
          if( idt.size() != size_t(L.nout) ) fail(what, ": inconsistent idt size");
          L.idt = ctx.weights.size(); ctx.weights.insert( ctx.weights.end(), idt.begin(), idt.end() );
        }
        if( !ctx.layers.empty() && int(ctx.layers.size()) > net.first_layer && ctx.layers.back().nout != L.nin ) fail(what, ": layer input does not match previous output");
        ctx.layers.push_back( L );
      }
      return net;
    }

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(

Reads a DeePMD model exported by exaStamp/scripts/python/deepmd_export.py (from any DeePMD-kit
model: .pth, .pt2, .dp, .yaml) and builds the tables used by deepmd_force, then raises rcut_max
to the model cutoff. exaStamp species are matched to the model's type_map by name.
Supported descriptors: se_e2_a, se_e2_r, se_e3, se_atten (DPA1, se_atten_v2; host only), se_t_tebd (se_e3_tebd), hybrid
(host only), dpa3 (repflows), dpa2 (repinit + three-body + repformers) and dpa4 (SeZM), message passing: host only, need a ghost update
before deepmd_force. Energies are in eV and lengths in Angstrom in the model file.

Place it in init_parameters, after species.

Usage example:

init_parameters:
  - species
  - deepmd_init:
      parameters: { file: "water_se_e2_a.dpnative.yaml" }

)EOF";
    }

    // DPA3 (repflows) descriptor of node dn; returns the descriptor size
    int read_dpa3(const YAML::Node& dn, DeepmdContext& ctx, int nt)
    {
      auto mp = std::make_shared<DPA3Model>();
      DPA3Model& m = *mp;
      m.nt = nt;
      m.n_dim = dn["n_dim"].as<int>(); m.e_dim = dn["e_dim"].as<int>(); m.a_dim = dn["a_dim"].as<int>();
      m.axis = dn["axis_neuron"].as<int>();
      m.nnei = dn["sel"].as< std::vector<int> >().at(0);
      m.a_sel = dn["a_sel"].as<int>();
      m.nca = dn["n_a_compress_dim"].as<int>(); m.eca = dn["e_a_compress_dim"].as<int>();
      m.compress = dn["a_compress_rate"].as<int>() == 0 ? 0 : ( dn["a_compress_use_split"].as<bool>() ? 1 : 2 );
      if( m.compress == 0 ) { m.nca = m.n_dim; m.eca = m.e_dim; }
      m.rc = dn["rcut"].as<double>(); m.rs = dn["rcut_smth"].as<double>();
      m.arc = dn["a_rcut"].as<double>(); m.ars = dn["a_rcut_smth"].as<double>();
      m.prot = dn["protection"].as<double>();
      m.exp_sw = dn["use_exp_switch"].as<bool>(); m.use_dist = dn["edge_init_use_dist"].as<bool>();
      m.smooth_edge = dn["smooth_edge_update"].as<bool>(); m.concat_tebd = dn["concat_output_tebd"].as<bool>();
      const std::string act = dn["activation"].as<std::string>();
      static const char* names[] = { "silut", "silu", "tanh", "relu", "gelu", "sigmoid", "softplus", "linear" };
      m.act = -1;
      for(int k = 0; k < 8; k++) if( act == names[k] ) m.act = k;
      if( act == "gelu_tf" ) m.act = DPA3_GELU;
      if( act == "none" ) m.act = DPA3_LINEAR;
      if( m.act < 0 ) fail("DPA3: unsupported activation '", act, "'");
      m.act_T = dn["act_threshold"].as<double>();
      { const double s = dpa3_sig( m.act_T ); m.act_slope = s + m.act_T * s * ( 1.0 - s ); m.act_const = m.act_T * s; }
      if( dn["use_dynamic_sel"].as<bool>() )
      {
        const double f = dn["sel_reduce_factor"].as<double>();
        m.s_sym = std::pow( m.nnei / f, -0.5 ); m.s_msg = 1.0 / ( m.nnei / f ); m.s_ang = std::pow( m.a_sel / f, -0.5 );
      }
      else { m.s_sym = std::pow( double(m.nnei), -0.5 ); m.s_msg = 1.0 / m.nnei; m.s_ang = std::pow( double(m.a_sel), -0.5 ); }
      auto put = [&](const std::vector<double>& v) { const size_t o = m.W.size(); m.W.insert( m.W.end(), v.begin(), v.end() ); return o; };
      auto lin = [&](const YAML::Node& x, int nin, int nout, const char* what) -> DPA3Lin
      {
        DPA3Lin L;
        if( !x || x.IsNull() ) return L;
        L.nin = x["nin"].as<int>(); L.nout = x["nout"].as<int>();
        if( ( nin >= 0 && L.nin != nin ) || ( nout >= 0 && L.nout != nout ) ) fail("DPA3 ", what, ": ", L.nin, " x ", L.nout, " instead of ", nin, " x ", nout);
        const auto w = x["w"].as< std::vector<double> >();
        if( w.size() != size_t(L.nin) * L.nout ) fail("DPA3 ", what, ": inconsistent weights");
        L.w = put( w );
        if( x["b"] && !x["b"].IsNull() ) { const auto b = x["b"].as< std::vector<double> >(); if( int(b.size()) != L.nout ) fail("DPA3 ", what, ": bias size"); L.b = put( b ); }
        return L;
      };
      const int nd = m.n_dim, ed = m.e_dim, ad = m.a_dim;
      m.edge_embd = lin( dn["edge_embd"], 1, ed, "edge_embd" );
      m.angle_embd = lin( dn["angle_embd"], 1, ad, "angle_embd" );
      const int adim = ad + m.nca + 2 * m.eca;
      for(const auto& x : dn["layers"])
      {
        DPA3Layer P;
        P.node_self = lin( x["node_self"], nd, nd, "node_self" );
        P.node_sym = lin( x["node_sym"], m.axis * ( nd + ed ), nd, "node_sym" );
        P.node_edge = lin( x["node_edge"], 2 * nd + ed, nd, "node_edge" );
        P.edge_self = lin( x["edge_self"], 2 * nd + ed, ed, "edge_self" );
        P.edge_angle1 = lin( x["edge_angle1"], adim, ed, "edge_angle1" );
        P.edge_angle2 = lin( x["edge_angle2"], ed, ed, "edge_angle2" );
        P.angle_self = lin( x["angle_self"], adim, ad, "angle_self" );
        if( m.compress == 2 )
        {
          P.compress_n = lin( x["compress_n"], nd, m.nca, "compress_n" );
          P.compress_e = lin( x["compress_e"], ed, m.eca, "compress_e" );
        }
        const auto rn = x["n_residual"].as< std::vector< std::vector<double> > >();
        const auto re = x["e_residual"].as< std::vector< std::vector<double> > >();
        const auto ra = x["a_residual"].as< std::vector< std::vector<double> > >();
        if( rn.size() != 3 || re.size() != 2 || ra.size() != 1 ) fail("DPA3: unexpected residual weights (res_residual, update_angle, one edge message)");
        for(int k = 0; k < 3; k++) P.n_res[k] = put( rn[k] );
        for(int k = 0; k < 2; k++) P.e_res[k] = put( re[k] );
        P.a_res = put( ra[0] );
        m.layers.push_back( P );
      }
      const auto tebd = dn["tebd"].as< std::vector< std::vector<double> > >();
      const auto davg = dn["davg"].as< std::vector< std::vector<double> > >();
      const auto dstd = dn["dstd"].as< std::vector< std::vector<double> > >();
      if( int(tebd.size()) != nt || int(davg.size()) != nt || int(dstd.size()) != nt ) fail("DPA3: tables for ", nt, " types expected");
      for(int t = 0; t < nt; t++)
      {
        if( int(tebd[t].size()) != nd ) fail("DPA3: type embedding width ", tebd[t].size(), " != n_dim");
        m.tebd.insert( m.tebd.end(), tebd[t].begin(), tebd[t].end() );
        m.davg.insert( m.davg.end(), davg[t].begin(), davg[t].end() );
        m.dstd.insert( m.dstd.end(), dstd[t].begin(), dstd[t].end() );
      }
      ctx.descriptor = DP_DPA3;
      ctx.ntypes = nt;
      ctx.rcut = m.rc;
      ctx.sel = { m.nnei };
      ctx.dim = nd + ( m.concat_tebd ? nd : 0 );
      ctx.dpa3 = mp;
      ctx.dpa3_buf = std::make_shared<DPA3Buffers>();
      return ctx.dim;
    }

    // one (non-hybrid) descriptor of node dn into ctx (its own weights and layers); returns the descriptor size
    int read_descriptor(const YAML::Node& dn, DeepmdContext& ctx, int nt)
    {
      ctx.ntypes = nt;
      const std::string desc = dn["descriptor"].as<std::string>();
      int descriptor = -1;
      if( desc == "se_e2_a" ) descriptor = DP_SE_E2_A;
      if( desc == "se_e2_r" ) descriptor = DP_SE_E2_R;
      if( desc == "se_e3" ) descriptor = DP_SE_E3;
      if( desc == "se_atten" ) descriptor = DP_SE_ATTEN;
      if( desc == "se_t_tebd" ) descriptor = DP_SE_T_TEBD;
      if( descriptor < 0 ) fail("unsupported descriptor '", desc, "' (se_e2_a, se_e2_r, se_e3, se_atten, se_t_tebd)");
      const bool atten = descriptor == DP_SE_ATTEN;
      const bool ttebd = descriptor == DP_SE_T_TEBD;
      const bool mixed = atten || ttebd;    // one neighbor list of the nnei nearest atoms, type embeddings

      ctx.descriptor = descriptor;
      ctx.ntypes = nt;
      ctx.rcut = dn["rcut"].as<double>();
      ctx.rcut_smth = dn["rcut_smth"].as<double>();
      ctx.protection = dn["protection"].as<double>();
      ctx.axis_neuron = dn["axis_neuron"].as<int>();
      ctx.type_one_side = dn["type_one_side"].as<bool>();

      const auto sel = dn["sel"].as< std::vector<int> >();
      if( int(sel.size()) != ( mixed ? 1 : nt ) ) fail("sel has ", sel.size(), " entries for ", nt, " types");
      ctx.mixed = mixed;
      ctx.sel.assign( sel.begin(), sel.end() );
      ctx.nnei = 0;
      for(int s : sel) ctx.nnei += s;

      const auto davg = dn["davg"].as< std::vector< std::vector< std::vector<double> > > >();
      const auto dstd = dn["dstd"].as< std::vector< std::vector< std::vector<double> > > >();
      for(int i = 0; i < nt; i++) for(int j = 0; j < nt; j++) for(int c = 0; c < 4; c++)
      {
        ctx.davg.push_back( davg.at(i).at(j).at(c) );
        ctx.dstd.push_back( dstd.at(i).at(j).at(c) );
      }


      const YAML::Node emb = dn["embedding"];
      if( descriptor == DP_SE_E2_R ) ctx.type_one_side = true;
      if( descriptor == DP_SE_E3 ) ctx.type_one_side = false;
      if( mixed )
      {
        const std::string mode = dn["tebd_input_mode"].as<std::string>();
        if( mode != "concat" && mode != "strip" ) fail("unsupported tebd_input_mode '", mode, "'");
        ctx.strip = mode == "strip";
        ctx.concat_output_tebd = dn["concat_output_tebd"].as<bool>();
        const auto tebd = dn["tebd"].as< std::vector< std::vector<double> > >();
        if( int(tebd.size()) != nt ) fail("tebd has ", tebd.size(), " rows for ", nt, " types");
        ctx.tebd_dim = tebd[0].size();
        for(const auto& row : tebd) { if( int(row.size()) != ctx.tebd_dim ) fail("inconsistent tebd"); ctx.tebd.insert( ctx.tebd.end(), row.begin(), row.end() ); }
      }
      if( ttebd ) { ctx.attn_smooth = dn["smooth"].as<bool>(); ctx.type_one_side = false; }
      if( atten )
      {
        const YAML::Node A = dn["attention"];
        ctx.attn_hidden = A["hidden"].as<int>();
        ctx.attn_scaling = A["scaling"].as<double>();
        ctx.attn_shift = A["shift"].as<double>();
        ctx.attn_normalize = A["normalize"].as<bool>();
        ctx.attn_smooth = A["smooth"].as<bool>();
        ctx.attn_dotr = A["dotr"].as<bool>();
        // the attention layers are read once ng (embedding output) is known, below
      }
      const size_t nemb = mixed ? ( ctx.strip ? 1 : ( ctx.type_one_side ? nt : nt * nt ) ) : ( ctx.type_one_side ? nt : nt * nt );
      if( emb.size() != nemb ) fail("expected ", nemb, " embedding networks, found ", emb.size());
      for(size_t k = 0; k < nemb; k++) ctx.embedding.push_back( read_net(emb[k], ctx, "embedding network " + std::to_string(k), DP_MAX_EMB_LAYERS, DP_MAX_EMB_WIDTH, "XSTAMP_DEEPMD_MAX_EMB") );

      ctx.ng = ctx.layers[ ctx.embedding[0].first_layer + ctx.embedding[0].nlayers - 1 ].nout;
      for(const auto& net : ctx.embedding)
      {
        if( ctx.layers[net.first_layer].nin != 1 ) fail("embedding network input must be 1");
        if( ctx.layers[net.first_layer + net.nlayers - 1].nout != ctx.ng ) fail("embedding networks have different output widths");
      }
      if( ctx.axis_neuron <= 0 || ctx.axis_neuron > ctx.ng ) fail("bad axis_neuron ", ctx.axis_neuron);
      if( ctx.ng * ctx.axis_neuron > DP_MAX_DESCRIPTOR ) fail("descriptor size ", ctx.ng * ctx.axis_neuron, " exceeds ", DP_MAX_DESCRIPTOR, " (XSTAMP_DEEPMD_MAX_DESCRIPTOR)");
      const int fit_nin = ctx.ng * ctx.axis_neuron + ( mixed && ctx.concat_output_tebd ? ctx.tebd_dim : 0 );
      if( ttebd && ctx.ng + ctx.tebd_dim > DP_MAX_EMB_WIDTH + 64 ) fail("se_t_tebd descriptor ng + tebd_dim exceeds ", DP_MAX_EMB_WIDTH + 64);
      if( atten )
      {
        const YAML::Node A = dn["attention"];
        const int ng = ctx.ng, hd = ctx.attn_hidden;
        auto put = [&](const YAML::Node& node, size_t expect, const char* what) -> size_t
        {
          const auto v = node.as< std::vector<double> >();
          if( v.size() != expect ) fail("attention ", what, ": ", v.size(), " values, ", expect, " expected");
          const size_t off = ctx.weights.size();
          ctx.weights.insert( ctx.weights.end(), v.begin(), v.end() );
          return off;
        };
        for(const auto& ln : A["layers"])
        {
          DPAttnLayer L;
          L.in_w = put( ln["in_w"], size_t(ng) * 3 * hd, "in_w" );
          L.in_b = put( ln["in_b"], size_t(3) * hd, "in_b" );
          L.out_w = put( ln["out_w"], size_t(hd) * ng, "out_w" );
          L.out_b = put( ln["out_b"], ng, "out_b" );
          L.ln_w = put( ln["ln_w"], ng, "ln_w" );
          L.ln_b = put( ln["ln_b"], ng, "ln_b" );
          L.ln_eps = ln["ln_eps"].as<double>();
          ctx.attn_layers.push_back( L );
        }
      }
      if( mixed && ctx.strip )
      {
        const int ng = ctx.ng;
        const auto T = dn["strip_table"].as< std::vector< std::vector<double> > >();
        if( T.size() != size_t( ctx.type_one_side ? nt : nt * nt ) ) fail("strip_table has ", T.size(), " rows");
        for(const auto& row : T) { if( int(row.size()) != ng ) fail("strip_table rows must have ", ng, " values"); ctx.strip_table.insert( ctx.strip_table.end(), row.begin(), row.end() ); }
      }

      // contribution of one empty neighbor slot to gr, Rh(pad) = -davg/dstd :
      // se_e2_a G(pad) x Rh(pad), se_e2_r G(pad) in column 0, se_e3 none (pairs with an empty slot have a zero angle term)
      if( descriptor == DP_SE_E3 || ttebd )
      {
        for(int i = 0; i < nt*nt; i++) for(int c = 1; c < 4; c++)
          if( ctx.davg[i*4 + c] != 0.0 ) fail("se_e3 with non-zero angular averages is not supported");
      }
      if( descriptor != DP_SE_E2_A && !atten && ctx.axis_neuron != 1 ) fail("axis_neuron must be 1 for ", desc);
      if( ttebd )
        for(int i = 0; i < nt; i++) for(int c = 1; c < 4; c++)
          if( ctx.davg[ (i*nt)*4 + c ] != 0.0 ) fail("se_t_tebd with non-zero angular averages is not supported");
      auto wk = std::make_unique<DPEmbeddingWork>();
      ctx.pad.assign( size_t(nt) * nt * ctx.ng * 4, 0.0 );
      // se_atten : one empty slot, row Rh_pad = -davg/dstd, type independent (zero averages in concat mode, or strip and
      // smooth: the type factor vanishes); masked in the attention, its feature goes through residual + out bias + norm
      for(int i = 0; i < nt && atten; i++)
      {
        double rh[4];
        bool zero = true;
        for(int c = 0; c < 4; c++) { rh[c] = - ctx.davg[ (i*nt)*4 + c ] / ctx.dstd[ (i*nt)*4 + c ]; zero = zero && rh[c] == 0.0; }
        if( zero ) continue;
        if( !ctx.strip || !ctx.attn_smooth ) fail("se_atten empty neighbor slots depend on the neighbor type (non-zero davg): not supported");
        const double* G = dp_mlp_forward( ctx.embedding[0], ctx.layers.data(), ctx.weights.data(), rh, *wk );
        std::vector<double> x( G, G + ctx.ng );
        for(const auto& L : ctx.attn_layers)
        {
          double mu = 0.0, var = 0.0;
          for(int c = 0; c < ctx.ng; c++) { x[c] += ctx.weights[L.out_b + c]; mu += x[c]; }
          mu /= ctx.ng;
          for(int c = 0; c < ctx.ng; c++) var += ( x[c] - mu ) * ( x[c] - mu );
          var /= ctx.ng;
          for(int c = 0; c < ctx.ng; c++) x[c] = ( x[c] - mu ) / std::sqrt( var + L.ln_eps ) * ctx.weights[L.ln_w + c] + ctx.weights[L.ln_b + c];
        }
        for(int a = 0; a < ctx.ng; a++) for(int c = 0; c < 4; c++) ctx.pad[ ( size_t(i*nt) * ctx.ng + a ) * 4 + c ] = x[a] * rh[c];
      }
      for(int i = 0; i < nt && descriptor != DP_SE_E3 && !mixed; i++) for(int j = 0; j < nt; j++)
      {
        double rh[4];
        for(int c = 0; c < 4; c++) rh[c] = - ctx.davg[ (i*nt + j)*4 + c ] / ctx.dstd[ (i*nt + j)*4 + c ];
        const DPNet& net = ctx.embedding[ ctx.type_one_side ? j : i*nt + j ];
        const double* G = dp_mlp_forward( net, ctx.layers.data(), ctx.weights.data(), rh, *wk );
        for(int a = 0; a < ctx.ng; a++) for(int c = 0; c < 4; c++)
          ctx.pad[ ( size_t(i*nt + j) * ctx.ng + a ) * 4 + c ] = ( descriptor == DP_SE_E2_R ) ? ( c == 0 ? G[a] : 0.0 ) : G[a] * rh[c];
      }

      ctx.dim = fit_nin;
      return fit_nin;
    }

    // DPA4 descriptor of node dn, with its GLU fitting network; returns the descriptor size
    int read_dpa4(const YAML::Node& dn, DeepmdContext& ctx, int nt)
    {
      auto mp = std::make_shared<DPA4Model>();
      DPA4Model& M = *mp;
      M.nt = nt;
      M.C = dn["channels"].as<int>(); M.lmax = dn["lmax"].as<int>();
      if( M.lmax < 1 || M.lmax > 3 ) fail("DPA4: lmax must be 1, 2 or 3");
      M.L = M.lmax + 1; M.D = M.L * M.L; M.Dm = M.L + 2 * M.lmax;
      M.F_n = dn["n_focus"] ? dn["n_focus"].as<int>() : 1; M.W = M.F_n * M.C;
      M.focus_tau = dn["focus_tau"] ? dn["focus_tau"].as<double>() : 1.0; M.focus_ls = dn["focus_ls"] ? dn["focus_ls"].as<double>() : 0.0;
      M.nnei = dn["sel"].as< std::vector<int> >().at(0);
      M.rc = dn["rcut"].as<double>(); M.eps = dn["eps"].as<double>(); M.weps = dn["wigner_eps"].as<double>(); M.deg_floor = dn["deg_floor"].as<double>();
      M.K = dn["frames"].as<int>(); M.Ng = dn["ng"].as<int>(); M.Ch = dn["ffn_channels"].as<int>(); M.Ce = M.K * M.Ch; M.rhid = dn["rhid"].as<int>();
      M.nk = dn["nk"].as<int>(); M.emb = dn["emb"].as<int>(); M.axis = dn["axis"].as<int>(); M.tdim = dn["tdim"].as<int>();
      M.hrf = dn["hrf"].as<int>(); M.s_r = dn["s_r"].as<int>(); M.s_gh = dn["s_gh"].as<int>(); M.fit_h = dn["fit_h"].as<int>();
      M.seed_eps = dn["seed_eps"].as<double>(); M.film_scale = dn["film_scale"].as<double>(); M.film_shift = dn["film_shift"].as<double>();
      M.fout_b = dn["fout_b"].as<double>();
      const int C = M.C, L = M.L, D = M.D, DK = D * M.K, K = M.K, Ch = M.Ch, Ce = M.Ce, lmax = M.lmax, LC = L * C, mC = lmax * C;
      auto arr = [&](const YAML::Node& x, size_t n, const char* what) { auto v = x.as< std::vector<double> >(); if( v.size() != n ) fail("DPA4 ", what, ": ", v.size(), " values, expected ", n); return v; };
      auto iarr = [&](const char* what, size_t n) { auto v = dn[what].as< std::vector<int> >(); if( n != size_t(-1) && v.size() != n ) fail("DPA4 ", what, ": ", v.size(), " values, expected ", n); return v; };
      M.env_c = dn["env_series"].as< std::vector<double> >(); M.rbf_c = dn["rbf_series"].as< std::vector<double> >();
      M.freq = dn["freq"].as< std::vector<double> >(); M.nrad = M.freq.size();
      const auto te = dn["te"].as< std::vector<double> >(), tt = dn["tt"].as< std::vector<double> >();
      if( te.size() < size_t(nt) * C || tt.size() < size_t(nt) * M.tdim ) fail("DPA4: type embeddings too small");
      M.te = te; M.tt = tt;
      M.rf1 = arr( dn["rf1"], size_t(M.nrad) * M.hrf, "rf1" ); M.rf2 = arr( dn["rf2"], size_t(M.hrf) * LC, "rf2" );
      if( dn["rf_scale"] ) { M.rf_norm.scale = arr( dn["rf_scale"], M.hrf, "rf_scale" ); M.rf_norm.eps = dn["rf_eps"].as<double>(); }
      const int Sc = M.s_r + 2 * M.tdim;
      M.s_r1 = arr( dn["s_r1"], size_t(M.nrad) * M.s_r, "s_r1" ); M.s_r2 = arr( dn["s_r2"], size_t(M.s_r) * M.s_r, "s_r2" );
      M.s_g1 = arr( dn["s_g1"], size_t(Sc) * M.s_gh, "s_g1" ); M.s_g2 = arr( dn["s_g2"], size_t(M.s_gh) * M.emb, "s_g2" );
      M.s_out = arr( dn["s_out"], size_t(M.emb) * M.axis * 2 * C, "s_out" );
      { const auto a = dn["s_act"].as< std::vector<bool> >(); for(int k = 0; k < 5; k++) M.s_a[k] = a.at(k); }
      if( dn["film_s"] ) { M.film_s.scale = arr( dn["film_s"], C, "film_s" ); M.film_s.eps = dn["film_s_eps"].as<double>(); }
      if( dn["film_f"] ) { M.film_f.scale = arr( dn["film_f"], C, "film_f" ); M.film_f.eps = dn["film_f_eps"].as<double>(); }
      M.C_l2 = arr( dn["C_l2"], 25 * 35, "C_l2" ); M.mono = iarr( "mono", 35 * 4 );
      M.nDq = M.lmax >= 3 ? 84 : 35;
      if( M.lmax >= 3 ) { M.C_l3 = arr( dn["C_l3"], 49 * 84, "C_l3" ); M.exp_l3 = iarr( "exp_l3", 84 * 4 ); }
      M.T = arr( dn["to_grid"], size_t(M.Ng) * DK, "to_grid" ); M.F = arr( dn["from_grid"], size_t(DK) * M.Ng, "from_grid" );
      M.rescale = arr( dn["rescale"], D, "rescale" ); M.sweight = arr( dn["sweight"], DK, "sweight" );
      M.coeff = iarr( "coeff", M.Dm ); M.deg_m = iarr( "deg_m", M.Dm ); M.sign_g = iarr( "gate_group", M.Dm );
      M.deg.resize( D ); for(int l = 0; l <= lmax; l++) for(int m = l*l; m < (l+1)*(l+1); m++) M.deg[m] = l;
      M.mix_c = iarr( "mix_compact", -1 ); M.mix_r = iarr( "mix_r", M.mix_c.size() ); M.mix_s = iarr( "mix_s", M.mix_c.size() );
      M.z_col = iarr( "z_offset", D - 1 ); M.z_slot = iarr( "z_slot", D - 1 );
      M.rot_r = iarr( "rot_r", -1 ); M.rot_d = iarr( "rot_d", M.rot_r.size() ); M.rot_o = iarr( "rot_o", M.rot_r.size() );
      auto eqn = [&](const YAML::Node& x, DPA4EqRms& N)
      {
        N.scale = arr( x["scale"], size_t(L) * C, "norm scale" ); N.bias = arr( x["bias"], C, "norm bias" ); N.bw = arr( x["bw"], D, "norm balance" ); N.eps = x["eps"].as<double>();
      };
      for(const auto& x : dn["blocks"])
      {
        DPA4Block P;
        const int F = M.F_n, W = M.W;
        // arrays with a focus axis: [rows][F][cols] -> [rows][cols] of focus f
        auto slice = [&](const std::vector<double>& a, size_t rows, size_t cols, int f)
        {
          std::vector<double> o( rows * cols );
          for(size_t r = 0; r < rows; r++) for(size_t c = 0; c < cols; c++) o[r * cols + c] = a[r * F * cols + f * cols + c];
          return o;
        };
        P.pre_mix = arr( x["pre_mix"], size_t(L) * C * W, "pre_mix" ); P.post_mix = arr( x["post_mix"], size_t(L) * W * C, "post_mix" );
        P.foc.resize( F );
        for(const auto& w : x["lin0"]) { const auto a = arr( w, size_t(LC) * F * LC, "lin0" ); for(int f = 0; f < F; f++) P.foc[f].lin0.push_back( slice( a, LC, LC, f ) ); }
        for(const auto& w : x["linm"]) { const auto a = arr( w, size_t(mC) * F * 2 * mC, "linm" ); for(int f = 0; f < F; f++) P.foc[f].linm.push_back( slice( a, mC, 2 * mC, f ) ); }
        for(const auto& w : x["gate"]) { const auto a = arr( w, size_t(C) * F * mC, "gate" ); for(int f = 0; f < F; f++) P.foc[f].gate.push_back( slice( a, C, mC, f ) ); }
        if( P.foc[0].gate.size() + 1 != P.foc[0].lin0.size() || P.foc[0].linm.size() != P.foc[0].lin0.size() ) fail("DPA4: inconsistent SO(2) layers");
        if( M.nk > 0 ) { P.mixw = arr( x["mixw"], size_t(L) * W * M.nk, "mixw" ); P.mixb = arr( x["mixb"], W, "mixb" ); }
        if( F > 1 ) P.rproj = arr( x["rproj"], size_t(C) * W, "rproj" );
        const auto qks = arr( x["qk_scale"], size_t(F) * C, "qk_scale" ), gts = arr( x["gate_scale"], size_t(F) * C, "gate_scale" );
        const auto qq = arr( x["q"], size_t(C) * F * C, "q" ), kk = arr( x["k"], size_t(C) * F * C, "k" );
        const auto lw = arr( x["logit_w"], size_t(C) * F, "logit_w" ), gw = arr( x["gate_w"], size_t(C) * F, "gate_w" );
        const auto ggt = arr( x["ggate"], size_t(2 * C) * F * C, "ggate" ), gr = arr( x["gres"], size_t(F) * C, "gres" );
        std::vector<double> lnull;
        if( x["log_null"].IsSequence() ) lnull = x["log_null"].as< std::vector<double> >(); else lnull.assign( 1, x["log_null"].as<double>() );
        if( int(lnull.size()) != F ) fail("DPA4: bad log_null");
        P.compete = F > 1 && x["comp_w"];
        for(int f = 0; f < F; f++)
        {
          DPA4Focus& Pf = P.foc[f];
          Pf.qk_norm.scale = slice( qks, 1, C, f ); Pf.qk_norm.eps = x["qk_eps"].as<double>();
          Pf.gate_norm.scale = slice( gts, 1, C, f ); Pf.gate_norm.eps = x["gate_eps"].as<double>();
          Pf.q = slice( qq, C, C, f ); Pf.k = slice( kk, C, C, f );
          Pf.logit_w = slice( lw, C, 1, f ); Pf.gate_w = slice( gw, C, 1, f ); Pf.log_null = lnull[f];
          Pf.ggate = slice( ggt, 2 * C, C, f ); Pf.gres = slice( gr, 1, C, f );
          if( P.compete )
          {
            Pf.comp_w = slice( arr( x["comp_w"], size_t(C) * F, "comp_w" ), C, 1, f );
            if( x["comp_scale"] ) { Pf.comp_norm.scale = slice( arr( x["comp_scale"], size_t(F) * C, "comp_scale" ), 1, C, f ); Pf.comp_norm.eps = x["comp_eps"].as<double>(); }
          }
        }
        P.gexp = arr( x["gexp"], size_t(L) * C * K * C, "gexp" ); P.gcon = arr( x["gcon"], size_t(L) * K * C * C, "gcon" );
        eqn( x["post_norm"], P.post_norm ); eqn( x["ffn_norm"], P.ffn_norm );
        P.flin1 = arr( x["flin1"], size_t(L) * C * 2 * Ce, "flin1" ); P.flin2 = arr( x["flin2"], size_t(L) * Ce * C, "flin2" );
        P.fl = arr( x["fl"], size_t(Ch) * Ch, "fl" ); P.fr = arr( x["fr"], size_t(Ch) * Ch, "fr" ); P.fo = arr( x["fo"], size_t(Ch) * Ch, "fo" );
        P.fgate = arr( x["fgate"], size_t(2 * Ch) * Ch, "fgate" );
        M.blocks.push_back( std::move(P) );
      }
      M.rlin1 = arr( dn["rlin1"], size_t(L) * C * 2 * Ce, "rlin1" );
      M.rl = arr( dn["rl"], size_t(2 * Ch) * M.rhid, "rl" ); M.rr = arr( dn["rr"], size_t(2 * Ch) * M.rhid, "rr" ); M.ro = arr( dn["ro"], size_t(M.rhid) * Ch, "ro" );
      M.rgate = arr( dn["rgate"], size_t(2 * Ch) * Ch, "rgate" ); M.rlin2 = arr( dn["rlin2"], size_t(Ch) * C, "rlin2" );
      M.fw = arr( dn["fw"], size_t(C) * 2 * M.fit_h, "fw" ); M.fb = arr( dn["fb"], 2 * M.fit_h, "fb" ); M.fout = arr( dn["fout"], M.fit_h, "fout" );
      ctx.dpa4 = mp;
      ctx.dpa4_buf = std::make_shared<DPA4Buffers>();
      ctx.descriptor = DP_DPA4;
      ctx.ntypes = nt;
      ctx.rcut = M.rc; ctx.nnei = M.nnei; ctx.sel = { M.nnei };
      return C;
    }

    // DPA2 descriptor of node dn (repinit and three-body sub-blocks read as se_atten / se_t_tebd contexts); returns the
    // descriptor size
    int read_dpa2(const YAML::Node& dn, DeepmdContext& ctx, int nt)
    {
      auto mp = std::make_shared<DPA2Model>();
      DPA2Model& M = *mp;
      DPA3Model& m = M.base;
      M.nt = nt; m.nt = nt;
      M.rep = std::make_shared<DeepmdContext>();
      M.dim_rep = read_descriptor( dn["repinit"], *M.rep, nt );
      if( M.rep->descriptor != DP_SE_ATTEN || !M.rep->attn_layers.empty() ) fail("DPA2: repinit must be an se_atten block without attention layers");
      if( dn["three_body"] && !dn["three_body"].IsNull() )
      {
        M.three = std::make_shared<DeepmdContext>();
        M.dim_3b = read_descriptor( dn["three_body"], *M.three, nt );
        if( M.three->descriptor != DP_SE_T_TEBD ) fail("DPA2: three-body block must be se_t_tebd");
      }
      M.dim0 = M.dim_rep + M.dim_3b;
      if( dn["repinit_dim"].as<int>() != M.dim0 ) fail("DPA2: repinit output ", M.dim0, " != ", dn["repinit_dim"].as<int>());
      const YAML::Node R = dn["repformer"];
      M.rc = R["rcut"].as<double>(); M.rs = R["rcut_smth"].as<double>(); M.prot = R["protection"].as<double>();
      M.nnei = R["sel"].as< std::vector<int> >().at(0);
      M.g1d = R["g1_dim"].as<int>(); M.g2d = R["g2_dim"].as<int>(); M.axis = R["axis_neuron"].as<int>();
      M.nh = R["attn2_nhead"].as<int>(); M.hid = R["attn2_hidden"].as<int>(); M.gate = R["attn2_has_gate"].as<bool>();
      M.shift = R["attnw_shift"].as<double>();
      M.s_sym = R["use_sqrt_nnei"].as<bool>() ? std::pow( double(M.nnei), -0.5 ) : 1.0 / M.nnei;
      M.out_conv = R["g1_out_conv"].as<bool>(); M.out_mlp = R["g1_out_mlp"].as<bool>();
      M.conv = R["update_g1_has_conv"].as<bool>(); M.grrg = R["update_g1_has_grrg"].as<bool>(); M.drrd = R["update_g1_has_drrd"].as<bool>();
      M.concat_tebd = dn["concat_output_tebd"].as<bool>();
      const std::string act = R["activation"].as<std::string>();
      static const char* names[] = { "silut", "silu", "tanh", "relu", "gelu", "sigmoid", "softplus", "linear" };
      m.act = -1;
      for(int k = 0; k < 8; k++) if( act == names[k] ) m.act = k;
      if( act == "gelu_tf" ) m.act = DPA3_GELU;
      if( act == "none" ) m.act = DPA3_LINEAR;
      if( m.act < 0 ) fail("DPA2: unsupported activation '", act, "'");
      if( R["act_threshold"] && !R["act_threshold"].IsNull() ) m.act_T = R["act_threshold"].as<double>();
      { const double s = dpa3_sig( m.act_T ); m.act_slope = s + m.act_T * s * ( 1.0 - s ); m.act_const = m.act_T * s; }
      auto put = [&](const std::vector<double>& v) { const size_t o = m.W.size(); m.W.insert( m.W.end(), v.begin(), v.end() ); return o; };
      auto lin = [&](const YAML::Node& x, int nin, int nout, const char* what, int bias) -> DPA3Lin
      {
        // bias: 0 none allowed, 1 required, -1 either
        DPA3Lin L;
        if( !x || x.IsNull() ) return L;
        L.nin = x["nin"].as<int>(); L.nout = x["nout"].as<int>();
        if( ( nin >= 0 && L.nin != nin ) || ( nout >= 0 && L.nout != nout ) ) fail("DPA2 ", what, ": ", L.nin, " x ", L.nout, " instead of ", nin, " x ", nout);
        const auto w = x["w"].as< std::vector<double> >();
        if( w.size() != size_t(L.nin) * L.nout ) fail("DPA2 ", what, ": inconsistent weights");
        L.w = put( w );
        if( x["b"] && !x["b"].IsNull() ) { const auto b = x["b"].as< std::vector<double> >(); if( int(b.size()) != L.nout ) fail("DPA2 ", what, ": bias size"); L.b = put( b ); }
        if( bias == 0 && L.b != size_t(-1) ) fail("DPA2 ", what, ": unexpected bias");
        if( bias == 1 && L.b == size_t(-1) ) fail("DPA2 ", what, ": bias expected");
        return L;
      };
      auto need = [&](const DPA3Lin& L, const char* what) { if( !L.valid() ) fail("DPA2: ", what, " missing"); };
      const int g1d = M.g1d, g2d = M.g2d, ax = M.axis, nh = M.nh, hid = M.hid;
      M.g2_embd = lin( R["g2_embd"], 1, g2d, "g2_embd", -1 ); need( M.g2_embd, "g2_embd" );
      M.shape = lin( dn["g1_shape"], M.dim0, g1d, "g1_shape", -1 );
      if( !M.shape.valid() && M.dim0 != g1d ) fail("DPA2: no g1 shape transform but repinit output ", M.dim0, " != g1_dim ", g1d);
      const auto tebd = dn["tebd"].as< std::vector< std::vector<double> > >();
      if( int(tebd.size()) != nt ) fail("DPA2: tebd for ", nt, " types expected");
      M.tdim = tebd[0].size();
      for(const auto& row : tebd) M.tebd.insert( M.tebd.end(), row.begin(), row.end() );
      M.tebd_tr = lin( dn["tebd_transform"], M.tdim, g1d, "tebd_transform", -1 );
      const auto davg = R["davg"].as< std::vector< std::vector<double> > >();
      const auto dstd = R["dstd"].as< std::vector< std::vector<double> > >();
      if( int(davg.size()) != nt || int(dstd.size()) != nt ) fail("DPA2: repformer statistics for ", nt, " types expected");
      for(int t = 0; t < nt; t++)
      {
        if( davg[t].size() != 4 || dstd[t].size() != 4 ) fail("DPA2: repformer statistics rows of 4 expected");
        for(int c = 1; c < 4; c++) if( davg[t][c] != 0.0 ) fail("DPA2: non-zero repformer angular averages are not supported");
        M.davg.insert( M.davg.end(), davg[t].begin(), davg[t].end() );
        M.dstd.insert( M.dstd.end(), dstd[t].begin(), dstd[t].end() );
      }
      const int dmlp = ( M.out_mlp ? 0 : g1d ) + ( M.conv && !M.out_conv ? g2d : 0 ) + ( M.grrg ? ax * g2d : 0 ) + ( M.drrd ? ax * g1d : 0 );
      const int nupd = ( M.out_mlp ? 1 : 0 ) + ( M.conv && M.out_conv ? 1 : 0 ) + 1;
      for(const auto& x : R["layers"])
      {
        DPA2Layer P;
        P.chnnl2 = x["update_chnnl_2"].as<bool>();
        P.attn = P.chnnl2 && x["update_g2_has_attn"].as<bool>();
        P.linear1 = lin( x["linear1"], dmlp, g1d, "linear1", 1 );
        if( M.out_mlp ) { P.self_mlp = lin( x["g1_self_mlp"], g1d, g1d, "g1_self_mlp", 1 ); need( P.self_mlp, "g1_self_mlp" ); }
        if( M.conv )
        {
          P.proj = M.out_conv ? lin( x["proj_g1g2"], g2d, g1d, "proj_g1g2", 0 ) : lin( x["proj_g1g2"], g1d, g2d, "proj_g1g2", 0 );
          need( P.proj, "proj_g1g2" );
        }
        if( P.chnnl2 ) { P.linear2 = lin( x["linear2"], g2d, g2d, "linear2", 1 ); need( P.linear2, "linear2" ); }
        if( P.attn )
        {
          P.mapqk = lin( x["mapqk"], g2d, hid * 2 * nh, "mapqk", 0 );
          P.mapv = lin( x["mapv"], g2d, g2d * nh, "mapv", 0 );
          P.head_map = lin( x["head_map"], g2d * nh, g2d, "head_map", 1 );
          need( P.mapqk, "mapqk" ); need( P.mapv, "mapv" ); need( P.head_map, "head_map" );
          const YAML::Node lm = x["attn2_lm"];
          const auto w = lm["w"].as< std::vector<double> >(), b = lm["b"].as< std::vector<double> >();
          if( int(w.size()) != g2d || int(b.size()) != g2d ) fail("DPA2: attn2_lm size");
          P.lnw = put( w ); P.lnb = put( b ); P.lneps = lm["eps"].as<double>();
        }
        const auto r1 = x["g1_residual"].as< std::vector< std::vector<double> > >();
        const auto r2 = x["g2_residual"].as< std::vector< std::vector<double> > >();
        if( int(r1.size()) != nupd ) fail("DPA2: ", r1.size(), " g1 residual weights, ", nupd, " expected");
        if( P.chnnl2 && int(r2.size()) != ( P.attn ? 2 : 1 ) ) fail("DPA2: unexpected g2 residual weights");
        for(const auto& r : r1) { if( int(r.size()) != g1d ) fail("DPA2: g1 residual size"); P.g1_res.push_back( put( r ) ); }
        for(const auto& r : r2) { if( int(r.size()) != g2d ) fail("DPA2: g2 residual size"); P.g2_res.push_back( put( r ) ); }
        M.layers.push_back( P );
      }
      if( M.layers.empty() ) fail("DPA2: no repformer layer");
      ctx.descriptor = DP_DPA2;
      ctx.ntypes = nt;
      ctx.rcut = std::max( { M.rc, M.rep->rcut, M.three ? M.three->rcut : 0.0 } );
      ctx.sel = { M.nnei };
      ctx.dim = g1d + ( M.concat_tebd ? M.tdim : 0 );
      ctx.dpa2 = mp;
      ctx.dpa2_buf = std::make_shared<DPA2Buffers>();
      return ctx.dim;
    }

    inline void execute() override final
    {
      m_fname = onika::data_file_path( parameters->file );
      YAML::Node doc;
      try { doc = YAML::LoadFile( m_fname ); }
      catch( const YAML::Exception& e ) { fail("cannot read model file: ", e.what()); }

      if( !doc["format"] || doc["format"].as<std::string>() != "exastamp-deepmd" ) fail("not a model exported by deepmd_export.py");
      if( doc["version"].as<int>() != 1 ) fail("unsupported model file version ", doc["version"].as<int>());

      DeepmdContext& ctx = *deepmd_ctx;
      ctx = DeepmdContext{};
      const auto type_map = doc["type_map"].as< std::vector<std::string> >();
      const int nt = type_map.size();
      if( nt == 0 || nt > DP_MAX_TYPES ) fail(nt, " types, at most ", DP_MAX_TYPES, " (XSTAMP_DEEPMD_MAX_TYPES)");
      const std::string desc = doc["descriptor"].as<std::string>();
      int dim = 0;
      if( desc == "hybrid" )
      {
        // sub-descriptors, each with its own neighbor selection within its rcut; descriptors concatenated
        ctx.descriptor = DP_HYBRID;
        ctx.ntypes = nt;
        ctx.rcut = 0.0;
        for(const auto& item : doc["list"])
        {
          auto sub = std::make_shared<DeepmdContext>();
          dim += read_descriptor( item, *sub, nt );
          ctx.rcut = std::max( ctx.rcut , sub->rcut );
          ctx.subs.push_back( sub );
        }
        if( ctx.subs.empty() ) fail("hybrid descriptor without sub-descriptors");
      }
      else if( desc == "dpa3" ) dim = read_dpa3( doc, ctx, nt );
      else if( desc == "dpa2" ) dim = read_dpa2( doc, ctx, nt );
      else if( desc == "dpa4" ) dim = read_dpa4( doc, ctx, nt );
      else dim = read_descriptor( doc, ctx, nt );
      ctx.dim = dim;

      // fitting networks, energy bias (top level)
      ctx.fitting_mixed_types = doc["fitting_mixed_types"].as<bool>();
      const auto bias = doc["energy_bias"].as< std::vector<double> >();
      if( int(bias.size()) != nt ) fail("energy_bias has ", bias.size(), " entries for ", nt, " types");
      ctx.energy_bias.assign( bias.begin(), bias.end() );
      const YAML::Node fit = doc["fitting"];
      const size_t nfit = desc == "dpa4" ? 0 : ( ctx.fitting_mixed_types ? 1 : nt );   // dpa4: GLU network in its model
      if( fit.size() != nfit ) fail("expected ", nfit, " fitting networks, found ", fit.size());
      for(size_t k = 0; k < nfit; k++) ctx.fitting.push_back( read_net(fit[k], ctx, "fitting network " + std::to_string(k), DP_MAX_FIT_LAYERS, DP_MAX_FIT_WIDTH, "XSTAMP_DEEPMD_MAX_FIT") );
      for(const auto& net : ctx.fitting)
      {
        if( ctx.layers[net.first_layer].nin != dim ) fail("fitting network input must be ", dim);
        if( ctx.layers[net.first_layer + net.nlayers - 1].nout != 1 ) fail("fitting network output must be 1");
      }
      if( desc == "hybrid" && dim > DP_MAX_DESCRIPTOR ) fail("hybrid descriptor size ", dim, " exceeds ", DP_MAX_DESCRIPTOR, " (XSTAMP_DEEPMD_MAX_DESCRIPTOR)");

      // exaStamp species -> model types, by name
      const auto& sp = *species;
      ctx.species_type.assign( sp.size(), -1 );
      for(size_t s = 0; s < sp.size(); s++)
      {
        const auto it = std::find( type_map.begin(), type_map.end(), std::string(sp[s].m_name) );
        if( it == type_map.end() )
        {
          std::string names;
          for(const auto& n : type_map) names += " " + n;
          fail("species '", sp[s].m_name, "' is not a type of the model (types:", names, ")");
        }
        ctx.species_type[s] = it - type_map.begin();
      }
      for(auto& sub : ctx.subs) sub->species_type = ctx.species_type;
      if( ctx.dpa4 ) ctx.dpa4->ebias.assign( ctx.energy_bias.begin(), ctx.energy_bias.end() );
      if( ctx.dpa2 ) { ctx.dpa2->rep->species_type = ctx.species_type; if( ctx.dpa2->three ) ctx.dpa2->three->species_type = ctx.species_type; }

      *rcut_max = std::max( *rcut_max , ctx.rcut );
      auto describe = [&](const DeepmdContext& c, const std::string& name)
      {
        lout << name << " rcut = " << c.rcut << " ang, sel =";
        for(int s : c.sel) lout << " " << s;
        lout << ", descriptor " << c.dim << ", " << c.weights.size() << " weights";
      };
      lout << "DeePMD " << desc << " : " << nt << " types";
      if( ctx.subs.empty() ) { lout << ", "; describe( ctx, "" ); }
      else for(size_t k = 0; k < ctx.subs.size(); k++) { lout << "\n  [" << k << "] "; describe( *ctx.subs[k], doc["list"][k]["descriptor"].as<std::string>() ); }
      lout << std::endl;
    }

  };

  ONIKA_AUTORUN_INIT(deepmd_init)
  {
    OperatorNodeFactory::instance()->register_factory("deepmd_init", make_simple_operator<DeepmdInit>);
  }

}
