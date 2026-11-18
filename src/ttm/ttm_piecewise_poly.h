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

#include <onika/cuda/cuda.h>
#include <onika/physics/units.h>
#include <onika/log.h>
#include <exaStamp/yaml_check_keys.h>
#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>

namespace exaStamp
{

  // Piecewise polynomial y(T) of the electronic temperature T [K], used for the Te-dependent electron
  // heat capacity Ce(Te) and electron-phonon coupling Gei(Te) of the TTM operator.
  //
  //   y(T) = unit_factor * y_scale * P_k( x_scale * T ),   P_k(x) = sum_i coeffs[k][i] * x^i
  //
  // piece k is used while T < breaks[k] (the last piece above the last break). Coefficients are given in
  // the fit's own unit (e.g. SI J/m^3/K), unit_factor converts to internal units once, on the host.
  //
  // Plain trivially-copyable struct with fixed-size arrays: captured by value in GPU functors.
  struct TtmPiecewisePoly
  {
    static constexpr int MaxPieces = 4;
    static constexpr int MaxCoeffs = 8;

    int    n_pieces = 0;                         // 0 => disabled (constant-parameter model)
    double breaks[MaxPieces-1] = {};             // [K], strictly increasing
    int    n_coeffs[MaxPieces] = {};
    double coeffs[MaxPieces][MaxCoeffs] = {};    // constant term first
    double x_scale = 1.0;
    double y_scale = 1.0;
    double unit_factor = 1.0;                    // fit unit -> internal units
    double integral_offset[MaxPieces] = {};      // integral of y from 0 to the lower bound of piece k (internal units * K)

    ONIKA_HOST_DEVICE_FUNC inline bool enabled() const { return n_pieces > 0; }

    ONIKA_HOST_DEVICE_FUNC inline int piece( double T ) const
    {
      int k = 0;
      while( k < n_pieces-1 && ! ( T < breaks[k] ) ) ++k;
      return k;
    }

    // lower bound of piece k (0 K for the first one: the integral starts at 0 K)
    ONIKA_HOST_DEVICE_FUNC inline double piece_lower_bound( int k ) const { return k==0 ? 0.0 : breaks[k-1]; }

    ONIKA_HOST_DEVICE_FUNC inline double poly( int k, double x ) const
    {
      double p = 0.0;
      for(int i=n_coeffs[k]-1;i>=0;i--) p = p*x + coeffs[k][i];
      return p;
    }

    // primitive of P_k, in x: sum_i c_i x^(i+1)/(i+1)
    ONIKA_HOST_DEVICE_FUNC inline double poly_primitive( int k, double x ) const
    {
      double p = 0.0;
      for(int i=n_coeffs[k]-1;i>=0;i--) p = p*x + coeffs[k][i] / double(i+1);
      return p * x;
    }

    // y(T) in internal units
    ONIKA_HOST_DEVICE_FUNC inline double value( double T ) const
    {
      return unit_factor * y_scale * poly( piece(T) , x_scale * T );
    }

    // integral of y from 0 K to T (internal units * K), e.g. the electronic energy density for Ce(Te).
    // Follows the pieces exactly, jumps at the breaks included.
    ONIKA_HOST_DEVICE_FUNC inline double integral( double T ) const
    {
      const int k = piece(T);
      const double scale = unit_factor * y_scale / x_scale;
      return integral_offset[k] + scale * ( poly_primitive( k , x_scale*T ) - poly_primitive( k , x_scale*piece_lower_bound(k) ) );
    }

    // Inverse of integral(): the T >= 0 such that integral(T) == E, assuming y > 0 on [0, T] (integral strictly
    // increasing; it is continuous across the breaks even where y jumps). Safeguarded Newton iteration
    // (d integral/dT = y) started from T_guess, falling back to bisection whenever a Newton step leaves the
    // current bracket [lo, hi] (jumps of y at the breaks, y ~ 0 near 0 K). E <= 0 gives 0.
    static constexpr int MaxInverseIterations = 100;
    ONIKA_HOST_DEVICE_FUNC inline double inverse_integral( double E, double T_guess ) const
    {
      if( ! ( E > 0.0 ) ) return 0.0;
      // Newton converges quadratically: after a Newton step s, the error is ~ (y'/2y) s^2 ~ s^2/T, so a
      // relative step below 1e-10 leaves T exact to roundoff. Bisection steps never end the iteration (their
      // error is of the order of the step), only a bracket down to a few ulps does.
      constexpr double newton_rtol = 1e-10;
      constexpr double bracket_rtol = 4e-16;
      double lo = 0.0;   // integral(lo) <= E
      double hi = -1.0;  // integral(hi) >= E; < 0 while no upper bound is known
      double T = T_guess > 0.0 ? T_guess : 1.0;
      for(int it=0; it<MaxInverseIterations; it++)
      {
        const double f = integral(T) - E;
        if( f == 0.0 ) return T;
        if( f < 0.0 ) lo = T; else hi = T;
        if( hi >= 0.0 && hi - lo <= bracket_rtol * hi ) return 0.5 * ( lo + hi );
        const double d = value(T);
        double Tn = ( d > 0.0 ) ? T - f/d : -1.0;
        bool newton = true;
        if( hi < 0.0 ) { if( ! ( Tn > lo ) ) { Tn = 2.0*lo + 1.0; newton = false; } } // no upper bound yet: grow
        else if( ! ( Tn > lo && Tn < hi ) ) { Tn = 0.5 * ( lo + hi ); newton = false; } // outside the bracket: bisect
        const double step = Tn > T ? Tn - T : T - Tn;
        if( newton && step <= newton_rtol * Tn ) return Tn;
        T = Tn;
      }
      return T;
    }

    // host only, after breaks/coeffs/scales are set
    inline void compute_integral_offsets()
    {
      const double scale = unit_factor * y_scale / x_scale;
      integral_offset[0] = 0.0;
      for(int k=1;k<n_pieces;k++)
      {
        integral_offset[k] = integral_offset[k-1]
                           + scale * ( poly_primitive( k-1 , x_scale*breaks[k-1] ) - poly_primitive( k-1 , x_scale*piece_lower_bound(k-1) ) );
      }
    }
  };

}

namespace YAML
{
  // Ce_te:
  //   unit: J/m^3/K            # unit of y_scale * P(x), converted to internal units
  //   x_scale: 1.0             # optional, x = x_scale * Te
  //   y_scale: 1.0             # optional
  //   breakpoints: [ 2200.0, 5500.0 ]                 # [K], one less than coefficients
  //   coefficients: [ [ 0.0, 67.6 ], [ ... ], [ ... ] ] # constant term first
  template<> struct convert<exaStamp::TtmPiecewisePoly>
  {
    static bool decode(const Node& node, exaStamp::TtmPiecewisePoly& v)
    {
      using P = exaStamp::TtmPiecewisePoly;
      if( ! node.IsMap() ) { return false; }
      exaStamp::check_yaml_keys( node , { "unit", "x_scale", "y_scale", "breakpoints", "coefficients" } , "TTM piecewise polynomial" );
      if( ! node["unit"] || ! node["coefficients"] )
      {
        onika::fatal_error() << "TTM piecewise polynomial: 'unit' and 'coefficients' are required" << std::endl;
      }

      v = P{};
      const auto coeffs = node["coefficients"].as< std::vector< std::vector<double> > >();
      const auto breaks = node["breakpoints"] ? node["breakpoints"].as< std::vector<double> >() : std::vector<double>{};
      if( coeffs.empty() || coeffs.size() > size_t(P::MaxPieces) )
      {
        onika::fatal_error() << "TTM piecewise polynomial: 1 to "<<P::MaxPieces<<" pieces expected, got "<<coeffs.size() << std::endl;
      }
      if( breaks.size()+1 != coeffs.size() )
      {
        onika::fatal_error() << "TTM piecewise polynomial: "<<coeffs.size()<<" pieces need "<<coeffs.size()-1<<" breakpoints, got "<<breaks.size() << std::endl;
      }
      v.n_pieces = int( coeffs.size() );
      for(size_t k=0;k<coeffs.size();k++)
      {
        if( coeffs[k].empty() || coeffs[k].size() > size_t(P::MaxCoeffs) )
        {
          onika::fatal_error() << "TTM piecewise polynomial: piece "<<k<<" needs 1 to "<<P::MaxCoeffs<<" coefficients, got "<<coeffs[k].size() << std::endl;
        }
        v.n_coeffs[k] = int( coeffs[k].size() );
        for(size_t i=0;i<coeffs[k].size();i++) v.coeffs[k][i] = coeffs[k][i];
      }
      for(size_t k=0;k<breaks.size();k++)
      {
        if( ( k==0 && ! ( breaks[0] > 0.0 ) ) || ( k>0 && ! ( breaks[k] > breaks[k-1] ) ) )
        {
          onika::fatal_error() << "TTM piecewise polynomial: breakpoints must be positive and strictly increasing" << std::endl;
        }
        v.breaks[k] = breaks[k];
      }
      if( node["x_scale"] ) v.x_scale = node["x_scale"].as<double>();
      if( node["y_scale"] ) v.y_scale = node["y_scale"].as<double>();
      if( ! ( v.x_scale > 0.0 ) )
      {
        onika::fatal_error() << "TTM piecewise polynomial: x_scale must be > 0" << std::endl;
      }
      v.unit_factor = onika::physics::make_quantity( 1.0 , node["unit"].as<std::string>() ).convert();
      v.compute_integral_offsets();
      return true;
    }
  };
}
