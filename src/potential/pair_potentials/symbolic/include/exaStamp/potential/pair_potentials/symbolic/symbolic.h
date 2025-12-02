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

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/cuda/cuda.h>
#include <exaStamp/unit_system.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/exprtk.hpp>
#include <cstring>

namespace exaStamp
{
  using namespace exanb;

  // Symbolic Parameters
  template<class ExpressionValueType>
  struct SymbolicParmsT
  {
    static inline constexpr size_t R_VARIABLE_INDEX = 0;
    static inline constexpr size_t MAX_VARIABLES = 8;
    using value_t = ExpressionValueType;
    using expression_t = exprtk::expression<value_t>;
    using symbol_table_t = exprtk::symbol_table<value_t>;
    
    struct ExpressionPrivData
    {
      symbol_table_t m_symbol_table;
      expression_t m_expression;
      value_t m_variables[MAX_VARIABLES];
      inline ExpressionPrivData() : m_symbol_table(symbol_table_t::e_mutable) {}
      ExpressionPrivData(const ExpressionPrivData&) =default;
      ExpressionPrivData(ExpressionPrivData &&) =default;
      ExpressionPrivData& operator = (const ExpressionPrivData&) =default;
      ExpressionPrivData& operator = (ExpressionPrivData &&) =default;
    };
    
    ExpressionPrivData * m_thread_ctx = nullptr;
    int m_max_threads = 0;
    
    inline void add_variable(const std::string& s, int var_index, const value_t& value, bool is_constant)
    {
      for(int i=0;i<m_max_threads;i++)
      {
        auto & ctx = m_thread_ctx[i];
        ctx.m_variables[var_index] = value;
        ctx.m_symbol_table.add_variable( s , ctx.m_variables[var_index] , is_constant );
      }
    }

    inline void compile(const std::string& expr)
    {
      for(int i=0;i<m_max_threads;i++)
      {
        auto & ctx = m_thread_ctx[i];
        exprtk::parser<value_t> parser;
        ctx.m_expression.register_symbol_table( ctx.m_symbol_table );
        parser.compile(expr,ctx.m_expression);
      }
    }
    
    inline void set_max_threads(int max_threads)
    {
      if( m_thread_ctx != nullptr ) delete [] m_thread_ctx;
      m_max_threads = max_threads;
      m_thread_ctx = new ExpressionPrivData[ m_max_threads ];
    }
    
  };
  
  using SymbolicParms = SymbolicParmsT<double>;

//# pragma omp declare simd uniform(p,ppp) notinbranch
  ONIKA_HOST_DEVICE_FUNC
  inline void symbolic_compute_energy(const SymbolicParms& p, const PairPotentialMinimalParameters&, double r, double& e, double& de)
  {
    assert( r > 0. );
    const int tid = omp_get_thread_num();
    assert( tid < p.m_max_threads );
    auto & ctx = p.m_thread_ctx[tid];
    ctx.m_variables[ p.R_VARIABLE_INDEX ] = r;
    e = ctx.m_expression.value();
    de = exprtk::derivative( ctx.m_expression ,"r");
  }
  
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::SymbolicParms;
  using onika::physics::Quantity;

  template<> struct convert<SymbolicParms>
  {
    static bool decode(const Node& node, SymbolicParms& v)
    {
      if( !node.IsMap() ) { return false; }
      if( !node["expression"] ) { return false; }
      if( !node["symbols"] ) { return false; }
      if( !node["symbols"].IsMap() ) { return false; }

      v.set_max_threads( omp_get_max_threads() );

      std::string expr_str = node["expression"].as<std::string>();
      
      int var_index = 0;
      v.add_variable( "r" , var_index, 0.0 , false );
      ++ var_index;
      
      for(auto sym : node["symbols"])
      {
        assert( var_index < SymbolicParms::MAX_VARIABLES );
        std::string sym_name = sym.first.as<std::string>();
        v.add_variable( sym_name, var_index , sym.second.as<Quantity>().convert() , true );
        ++ var_index;
      }
      v.compile(expr_str);

      return true;
    }
  };
}

