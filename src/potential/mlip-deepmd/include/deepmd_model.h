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

#pragma once

// Native DeePMD model: weights of the embedding and fitting networks, environment
// matrix statistics, and the dense layers (forward and hand-written backward).
// The model file is written by exaStamp/scripts/python/deepmd_export.py from any
// DeePMD-kit model; the reference implementation is deepmd/dpmodel (numpy).

#include <yaml-cpp/yaml.h>
#include <exaStamp/yaml_check_keys.h>
#include <onika/memory/allocator.h>
#include <onika/cuda/cuda.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <memory>
#include <vector>

// Size limits of the per-atom scratch. On GPU the scratch lives in block shared memory
// (see deepmd_force_op.h), so these bound the models a build can run. CMake cache variables.
#ifndef XSTAMP_DEEPMD_MAX_TYPES
#define XSTAMP_DEEPMD_MAX_TYPES 16
#endif
#ifndef XSTAMP_DEEPMD_MAX_EMB_WIDTH
#define XSTAMP_DEEPMD_MAX_EMB_WIDTH 128
#endif
#ifndef XSTAMP_DEEPMD_MAX_EMB_LAYERS
#define XSTAMP_DEEPMD_MAX_EMB_LAYERS 4
#endif
#ifndef XSTAMP_DEEPMD_MAX_FIT_WIDTH
#define XSTAMP_DEEPMD_MAX_FIT_WIDTH 256
#endif
#ifndef XSTAMP_DEEPMD_MAX_FIT_LAYERS
#define XSTAMP_DEEPMD_MAX_FIT_LAYERS 4
#endif
#ifndef XSTAMP_DEEPMD_MAX_DESCRIPTOR
#define XSTAMP_DEEPMD_MAX_DESCRIPTOR 2048
#endif

namespace exaStamp
{

  static constexpr int DP_MAX_TYPES = XSTAMP_DEEPMD_MAX_TYPES;
  static constexpr int DP_MAX_EMB_WIDTH = XSTAMP_DEEPMD_MAX_EMB_WIDTH;     // embedding layer width, embedding output (ng)
  static constexpr int DP_MAX_EMB_LAYERS = XSTAMP_DEEPMD_MAX_EMB_LAYERS;
  static constexpr int DP_MAX_FIT_WIDTH = XSTAMP_DEEPMD_MAX_FIT_WIDTH;     // fitting hidden layer width
  static constexpr int DP_MAX_FIT_LAYERS = XSTAMP_DEEPMD_MAX_FIT_LAYERS;   // fitting layers, output layer included
  static constexpr int DP_MAX_DESCRIPTOR = XSTAMP_DEEPMD_MAX_DESCRIPTOR;   // descriptor size ng x axis_neuron

  // User parameters (YAML) of deepmd_init
  struct DeepmdParams
  {
    std::string file; // model exported by deepmd_export.py
  };

  // descriptors (deepmd/dpmodel/descriptor)
  enum DPDescriptor : int
  {
    DP_SE_E2_A = 0,  // se_e2_a : two-body embedding of the full environment matrix row
    DP_SE_E2_R,      // se_e2_r : two-body, radial only
    DP_SE_E3,        // se_e3   : three-body, angles between neighbor pairs
    DP_SE_ATTEN,     // se_atten (DPA1, se_atten_v2) : mixed neighbor list, type embeddings, gated self-attention
    DP_SE_T_TEBD,    // se_e3_tebd / se_t_tebd : three-body angles over a mixed neighbor list, type embeddings
    DP_HYBRID,       // hybrid : concatenation of sub-descriptors (DeepmdContext::subs)
    DP_DPA3,         // dpa3 (repflows) : message passing between atoms, host (DeepmdContext::dpa3)
    DP_DPA4,         // dpa4 (SeZM) : SO(3) message passing between atoms, host (DeepmdContext::dpa4)
    DP_DPA2          // dpa2 (repinit + three-body + repformers) : message passing between atoms, host (DeepmdContext::dpa2)
  };

  struct DPA3Model;
  struct DPA3Buffers;
  struct DPA2Model;
  struct DPA2Buffers;
  struct DPA4Model;
  struct DPA4Buffers;

  enum DPActivation : int
  {
    DP_ACT_NONE = 0,
    DP_ACT_TANH,
    DP_ACT_RELU,
    DP_ACT_GELU,      // tanh approximation (deepmd 'gelu' / 'gelu_tf')
    DP_ACT_SOFTPLUS,
    DP_ACT_SIGMOID,
    DP_ACT_SILU
  };

  // y = act(x W + b) [* idt] [+ x | + (x,x)], offsets index DeepmdContext::weights
  struct DPLayer
  {
    int nin = 0;
    int nout = 0;
    int activation = DP_ACT_NONE;
    int resnet = 0;          // 0: none, 1: + x (nout == nin), 2: + (x,x) (nout == 2 nin)
    size_t w = 0;            // nin x nout, row-major
    size_t b = 0;            // nout
    size_t idt = size_t(-1); // nout, or none
  };

  struct DPNet
  {
    int first_layer = 0;
    int nlayers = 0;
  };

  // one se_atten attention layer: input projection ng -> 3 hidden (q, k, v), output projection hidden -> ng,
  // layer norm over ng (offsets in DeepmdContext::weights)
  struct DPAttnLayer
  {
    size_t in_w = 0, in_b = 0, out_w = 0, out_b = 0, ln_w = 0, ln_b = 0;
    double ln_eps = 1e-5;
  };

  // Output of deepmd_init, consumed by deepmd_force.
  // Model types are the entries of the model's type_map; exaStamp species map to them by name.
  struct DeepmdContext
  {
    onika::memory::CudaMMVector<double> weights;
    onika::memory::CudaMMVector<DPLayer> layers;
    onika::memory::CudaMMVector<DPNet> embedding;     // by neighbor type (one side), else type pair a x ntypes + b
    onika::memory::CudaMMVector<DPNet> fitting;       // by center type, or a single network
    onika::memory::CudaMMVector<double> davg;         // [center][neighbor][4]
    onika::memory::CudaMMVector<double> dstd;         // [center][neighbor][4]
    onika::memory::CudaMMVector<double> pad;          // [center][neighbor][ng][4] : contribution of one empty slot to gr
    onika::memory::CudaMMVector<double> energy_bias;  // [type], eV
    onika::memory::CudaMMVector<int> sel;             // [type]
    onika::memory::CudaMMVector<int> species_type;    // exaStamp species -> model type
    int descriptor = DP_SE_E2_A;
    int ntypes = 0;
    int ng = 0;          // embedding output width
    int axis_neuron = 0;
    int nnei = 0;        // sum of sel
    bool type_one_side = true;
    bool fitting_mixed_types = false;
    double rcut = 0.0;
    double rcut_smth = 0.0;
    double protection = 0.0;
    // se_atten
    onika::memory::CudaMMVector<DPAttnLayer> attn_layers;
    onika::memory::CudaMMVector<double> tebd;         // [type][tebd_dim]
    onika::memory::CudaMMVector<double> strip_table;  // strip mode: [neighbor type (one side) or center x ntypes + neighbor][ng]
    int tebd_dim = 0;
    int attn_hidden = 0;
    double attn_scaling = 1.0;
    double attn_shift = 20.0;
    bool attn_normalize = true;
    bool attn_smooth = true;
    bool attn_dotr = true;
    bool strip = false;
    bool concat_output_tebd = true;
    bool mixed = false;              // a single neighbor list of the nnei nearest atoms (se_atten, se_t_tebd)
    int dim = 0;                     // descriptor size (fitting network input)
    std::vector< std::shared_ptr<DeepmdContext> > subs;   // hybrid: sub-descriptors (no fitting networks)
    std::shared_ptr<DPA3Model> dpa3;                     // dpa3 model
    std::shared_ptr<DPA3Buffers> dpa3_buf;               // dpa3 per-run buffers, kept between runs
    std::shared_ptr<DPA2Model> dpa2;                     // dpa2 model (repinit / three-body contexts inside)
    std::shared_ptr<DPA2Buffers> dpa2_buf;               // dpa2 per-run buffers
    std::shared_ptr<DPA4Model> dpa4;                     // dpa4 model (with its own fitting network)
    std::shared_ptr<DPA4Buffers> dpa4_buf;               // dpa4 per-run buffers
  };

  // ------------------------------------------------------------------ activations

  ONIKA_HOST_DEVICE_FUNC inline double dp_sigmoid(double z) { return 1.0 / ( 1.0 + std::exp(-z) ); }

  ONIKA_HOST_DEVICE_FUNC inline double dp_activation(int a, double z)
  {
    switch( a )
    {
      case DP_ACT_TANH : return std::tanh(z);
      case DP_ACT_RELU : return z > 0.0 ? z : 0.0;
      case DP_ACT_GELU : { const double k = 0.7978845608028654; return 0.5 * z * ( 1.0 + std::tanh( k * ( z + 0.044715 * z*z*z ) ) ); }
      case DP_ACT_SOFTPLUS : return z > 30.0 ? z : std::log1p( std::exp(z) );
      case DP_ACT_SIGMOID : return dp_sigmoid(z);
      case DP_ACT_SILU : return z * dp_sigmoid(z);
      default : return z;
    }
  }

  ONIKA_HOST_DEVICE_FUNC inline double dp_activation_grad(int a, double z)
  {
    switch( a )
    {
      case DP_ACT_TANH : { const double t = std::tanh(z); return 1.0 - t*t; }
      case DP_ACT_RELU : return z > 0.0 ? 1.0 : 0.0;
      case DP_ACT_GELU :
      {
        const double k = 0.7978845608028654;
        const double t = std::tanh( k * ( z + 0.044715 * z*z*z ) );
        return 0.5 * ( 1.0 + t ) + 0.5 * z * ( 1.0 - t*t ) * k * ( 1.0 + 3.0 * 0.044715 * z*z );
      }
      case DP_ACT_SOFTPLUS : return dp_sigmoid(z);
      case DP_ACT_SIGMOID : { const double s = dp_sigmoid(z); return s * ( 1.0 - s ); }
      case DP_ACT_SILU : { const double s = dp_sigmoid(z); return s + z * s * ( 1.0 - s ); }
      default : return 1.0;
    }
  }

  // ------------------------------------------------------------------ dense networks

  // Activations kept by the forward pass for the backward pass.
  // x[l] is the output of layer l (the network input is passed separately), z[l] the pre-activation of layer l.
  template<int MaxLayers, int MaxWidth>
  struct DPMlpWork
  {
    static constexpr int max_layers = MaxLayers;
    static constexpr int max_width = MaxWidth;
    double x[MaxLayers][MaxWidth];
    double z[MaxLayers][MaxWidth];
  };

  // one thread per embedding network evaluation
  using DPEmbeddingWork = DPMlpWork<DP_MAX_EMB_LAYERS, DP_MAX_EMB_WIDTH>;

  // forward pass of net on input in[0..nin), returns a pointer to the output
  template<class WorkT>
  ONIKA_HOST_DEVICE_FUNC inline const double* dp_mlp_forward(const DPNet& net, const DPLayer* layers, const double* W,
                                                             const double* in, WorkT& wk)
  {
    const double* x = in;
    for(int l = 0; l < net.nlayers; l++)
    {
      const DPLayer& L = layers[net.first_layer + l];
      const double* w = W + L.w;
      const double* b = W + L.b;
      double* z = wk.z[l];
      double* y = wk.x[l];
      for(int o = 0; o < L.nout; o++) z[o] = b[o];
      for(int i = 0; i < L.nin; i++)
      {
        const double xi = x[i];
        const double* wi = w + size_t(i) * L.nout;
        for(int o = 0; o < L.nout; o++) z[o] += xi * wi[o];
      }
      for(int o = 0; o < L.nout; o++)
      {
        double v = dp_activation(L.activation, z[o]);
        if( L.idt != size_t(-1) ) v *= W[L.idt + o];
        if( L.resnet == 1 ) v += x[o];
        else if( L.resnet == 2 ) v += x[ o % L.nin ];
        y[o] = v;
      }
      x = y;
    }
    return x;
  }

  // backward pass (vector-Jacobian product): dout[0..nout) -> din[0..nin), after dp_mlp_forward on the same input.
  // g0 and g1 are scratch arrays as wide as the widest layer; dout may alias one of them.
  template<class WorkT>
  ONIKA_HOST_DEVICE_FUNC inline void dp_mlp_backward(const DPNet& net, const DPLayer* layers, const double* W, const double* in,
                                                     const WorkT& wk, const double* dout, double* din, double* g0, double* g1)
  {
    const double* dy = dout;
    for(int l = net.nlayers - 1; l >= 0; l--)
    {
      const DPLayer& L = layers[net.first_layer + l];
      const double* w = W + L.w;
      const double* z = wk.z[l];
      double* dz = ( dy == g0 ) ? g1 : g0;
      for(int o = 0; o < L.nout; o++)
      {
        double g = dy[o] * dp_activation_grad(L.activation, z[o]);
        if( L.idt != size_t(-1) ) g *= W[L.idt + o];
        dz[o] = g;
      }
      double* dx = ( l == 0 ) ? din : ( ( dz == g0 ) ? g1 : g0 );
      for(int i = 0; i < L.nin; i++)
      {
        const double* wi = w + size_t(i) * L.nout;
        double s = 0.0;
        for(int o = 0; o < L.nout; o++) s += wi[o] * dz[o];
        if( L.resnet == 1 ) s += dy[i];
        else if( L.resnet == 2 ) s += dy[i] + dy[i + L.nin];
        dx[i] = s;
      }
      dy = dx;
    }
    (void) in;
  }

}

namespace YAML
{
  template<> struct convert<exaStamp::DeepmdParams>
  {
    static bool decode(const Node& node, exaStamp::DeepmdParams& v)
    {
      if( !node.IsMap() ) { return false; }
      exaStamp::check_yaml_keys( node , { "file" } , "deepmd parameters" );
      if( !node["file"] ) { return false; }
      v.file = node["file"].as<std::string>();
      return true;
    }
  };
}
