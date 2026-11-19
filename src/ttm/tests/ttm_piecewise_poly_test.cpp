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

// Unit test of the Te-dependent TTM parameters (ttm_piecewise_poly.h) on the gold Ce(Te)/Gei(Te) fit:
// spot values, continuity at the switch points, positivity, array vs scalar evaluation, integral of Ce.
// The fits are parsed through the real YAML decoder; values are compared in SI (divided by unit_factor).

#include "ttm_piecewise_poly.h"

#include <cmath>
#include <cstdio>
#include <vector>

using exaStamp::TtmPiecewisePoly;

static int g_failures = 0;

static void check( bool ok, const char* what, double got, double expected )
{
  std::printf( "%-4s %-48s got=%.10g expected=%.10g\n", ok ? "OK" : "FAIL", what, got, expected );
  if( ! ok ) ++g_failures;
}

static void check_rel( const char* what, double got, double expected, double rtol )
{
  check( std::fabs(got-expected) <= rtol*std::fabs(expected), what, got, expected );
}

static const char* ce_yaml = R"EOF(
unit: J/m^3/K
breakpoints: [ 2200.0, 5500.0 ]
coefficients:
  - [ 0.0, 67.6 ]
  - [ 22135.0, 12.685, 1.8962e-2 ]
  - [ -4.027e5, 210.74, -3.0463e-3, 4.0076e-8, 0.0, 1.0652e-18 ]
)EOF";

static const char* gei_yaml = R"EOF(
unit: J/s/m^3/K
x_scale: 1.0e-4
y_scale: 1.0e17
breakpoints: [ 2100.0, 4400.0 ]
coefficients:
  - [ 0.251 ]
  - [ 0.25353, 0.042859, -0.92027, 3.1317 ]
  - [ -0.53893, 2.3867, -0.83833, 0.15177, -0.013335, -4.5018e-4 ]
)EOF";

int main()
{
  // exaStamp internal unit system, as set at run time from data/config/config_defaults.msp
  onika::physics::set_internal_unit_system( YAML::Load( "{ length: angstrom, mass: Dalton, time: picosecond, charge: elementary_charge, "
    "temperature: kelvin, amount: particle, luminosity: candela, angle: radian, energy: joule }" ).as<onika::physics::UnitSystem>() );

  const TtmPiecewisePoly ce_p = YAML::Load(ce_yaml).as<TtmPiecewisePoly>();
  const TtmPiecewisePoly gei_p = YAML::Load(gei_yaml).as<TtmPiecewisePoly>();
  auto Ce  = [&](double T) { return ce_p.value(T) / ce_p.unit_factor; };   // SI
  auto Gei = [&](double T) { return gei_p.value(T) / gei_p.unit_factor; }; // SI
  std::printf( "unit factors: Ce %.10g, Gei %.10g\n", ce_p.unit_factor, gei_p.unit_factor );
  // 1 J = 1/(Da*1e4) Da.ang^2/ps^2, 1/m^3 = 1e-30 /ang^3, 1/s = 1e-12 /ps
  const double J_internal = 1.0 / ( 1.66053906660e-27 * 1.0e4 );
  check_rel( "J/m^3/K   -> Da/(ang.ps^2.K)", ce_p.unit_factor, J_internal * 1e-30, 1e-8 ); // onika uses an older CODATA Dalton
  check_rel( "J/s/m^3/K -> Da/(ang.ps^3.K)", gei_p.unit_factor, J_internal * 1e-42, 1e-8 );

  // 1. spot values
  check_rel( "Ce(1000)", Ce(1000.0), 67600.0, 1e-12 );
  check_rel( "Ce(3000)", Ce(3000.0), 230848.0, 1e-12 );
  check_rel( "Gei(1000)", Gei(1000.0), 2.51e16, 1e-12 );
  check_rel( "Gei(3000)/1e17", Gei(3000.0)/1e17, 0.268119, 2e-6 );
  // high branches, against an independent reference implementation of the fits (numpy, double precision)
  check_rel( "Ce(8000)  vs python ref", Ce(8000.0), 1108810.6164736, 1e-13 );
  check_rel( "Ce(1e4)   vs python ref", Ce(1e4), 1440252.52, 1e-13 );
  check_rel( "Gei(8000) vs python ref", Gei(8000.0), 9.059955090175998e+16, 1e-13 );
  check_rel( "Gei(1e4)  vs python ref", Gei(1e4), 1.1474248199999995e+17, 1e-13 );

  // 2. continuity at the switch points (left/right values from the reference table)
  const double eps = 1e-6;
  struct Sw { const char* name; bool is_ce; double T, left, right, rtol; };
  const Sw sw[] = {
    { "Ce   2200 K", true , 2200.0, 148720.0, 141818.0, 1e-5 },
    { "Ce   5500 K", true , 5500.0, 665503.0, 670892.0, 1e-5 },
    { "Gei  2100 K", false, 2100.0, 0.25100 , 0.25095 , 1e-4 },
    { "Gei  4400 K", false, 4400.0, 0.36099 , 0.36134 , 1e-4 } };
  for(const auto& s : sw)
  {
    const double scale = s.is_ce ? 1.0 : 1e-17;
    const double l = ( s.is_ce ? Ce(s.T-eps) : Gei(s.T-eps) ) * scale;
    const double r = ( s.is_ce ? Ce(s.T+eps) : Gei(s.T+eps) ) * scale;
    char buf[64];
    std::snprintf(buf,sizeof(buf),"%s left",s.name);  check_rel( buf, l, s.left , s.rtol );
    std::snprintf(buf,sizeof(buf),"%s right",s.name); check_rel( buf, r, s.right, s.rtol );
    std::snprintf(buf,sizeof(buf),"%s jump [%%]",s.name);
    const double jump = 100.0*(r-l)/l;
    const double max_jump = s.is_ce ? ( s.T < 3000.0 ? 4.7 : 0.9 ) : 0.2;
    check( std::fabs(jump) <= max_jump, buf, jump, max_jump );
  }
  check( ce_p.piece(2200.0) == 1 && ce_p.piece(5500.0) == 2 && gei_p.piece(2100.0) == 1, "switch point belongs to the upper piece", 0, 0 );

  // 3. positivity on (0, 1e4] K
  bool positive = true; double worst = 1e300;
  for(int i=1;i<=100000;i++)
  {
    const double T = 1e4 * i / 100000.0;
    const double c = Ce(T), g = Gei(T);
    positive = positive && c > 0.0 && g > 0.0;
    worst = std::fmin( worst, std::fmin( c/Ce(1e4), g/Gei(1e4) ) );
  }
  check( positive, "Ce>0 and Gei>0 on (0,1e4] K", worst, 0.0 );

  // 4. array evaluation == element-wise scalar evaluation (bit-for-bit)
  std::vector<double> T(20001), ce_arr(T.size()), gei_arr(T.size());
  for(size_t i=0;i<T.size();i++) T[i] = 0.5 * double(i);
  const TtmPiecewisePoly ce_copy = ce_p, gei_copy = gei_p; // by-value copy, as captured by GPU functors
  for(size_t i=0;i<T.size();i++) { ce_arr[i] = ce_copy.value(T[i]); gei_arr[i] = gei_copy.value(T[i]); }
  bool same = true;
  for(size_t i=0;i<T.size();i++) same = same && ce_arr[i] == ce_p.value(T[i]) && gei_arr[i] == gei_p.value(T[i]);
  check( same, "array evaluation == scalar evaluation", 0, 0 );

  // 5. integral of Ce (electronic energy density) vs composite Simpson, across both jumps
  for(double Tmax : { 1000.0, 2200.0, 3000.0, 5500.0, 8000.0, 1e4 })
  {
    // integrate piece by piece so Simpson never straddles a jump
    const double bounds[] = { 0.0, 2200.0, 5500.0, 1e300 };
    double ref = 0.0;
    for(int k=0;k<3 && bounds[k]<Tmax;k++)
    {
      const double a = bounds[k], b = std::fmin(bounds[k+1],Tmax);
      const int n = 2000; const double h = (b-a)/n;
      // evaluate strictly inside [a,b): the piece is selected with '<', b itself belongs to the next piece
      auto f = [&](double t) { return ce_p.unit_factor*ce_p.y_scale*ce_p.poly(k,t); };
      double sacc = f(a) + f(b);
      for(int i=1;i<n;i++) sacc += ( i%2 ? 4.0 : 2.0 ) * f(a+i*h);
      ref += sacc*h/3.0;
    }
    char buf[64]; std::snprintf(buf,sizeof(buf),"integral_0^%g Ce dT (internal)",Tmax);
    check_rel( buf, ce_p.integral(Tmax), ref, 1e-12 );
  }
  check_rel( "integral_0^1000 Ce dT (SI) = 67.6/2*1e6", ce_p.integral(1000.0)/ce_p.unit_factor, 33.8e6, 1e-12 );

  // 6. inverse of the integral (energy-form Te update): Te such that integral_0^Te Ce dT = E
  //    a) linear Ce = gamma*T (first gold piece): E = gamma*T^2/2, T = sqrt(2E/gamma), analytic
  {
    const double gamma = 67.6 * ce_p.unit_factor;
    double worst = 0.0;
    for(double Tref : { 1.0, 10.0, 300.0, 1000.0, 2000.0, 2199.0 })
    {
      const double E = 0.5 * gamma * Tref * Tref;
      for(double g : { 0.0, 0.5*Tref, Tref*(1.0+1e-3), 3.0*Tref, 1e5 }) // bad guesses included
      {
        const double T = ce_p.inverse_integral( E, g );
        worst = std::fmax( worst, std::fabs( T - std::sqrt( 2.0*E/gamma ) ) / Tref );
      }
    }
    check( worst <= 4e-16, "inverse_integral, Ce = gamma*T, vs sqrt(2E/gamma)", worst, 4e-16 );
  }
  //    b) round trip T -> E -> T over the whole fit, both jumps included, with the explicit-step-like guess
  {
    double worst = 0.0, worst_T = 0.0;
    for(int i=1;i<=100000;i++)
    {
      const double T0 = 1e4 * i / 100000.0;
      for(double g : { T0, T0*(1.0+1e-2), T0*(1.0-1e-2), 0.0 })
      {
        const double T = ce_p.inverse_integral( ce_p.integral(T0), g );
        const double r = std::fabs( T - T0 ) / T0;
        if( r > worst ) { worst = r; worst_T = T0; }
      }
    }
    std::printf( "     worst round trip at T=%g K\n", worst_T );
    check( worst <= 1e-15, "inverse_integral(integral(T)) == T on (0,1e4] K", worst, 1e-15 );
    // around the breaks, where Ce jumps: E is continuous there, the inverse must be too
    double worst_b = 0.0;
    for(double Tb : { 2200.0, 5500.0 })
    for(double d : { -1e-6, -1e-9, 0.0, 1e-9, 1e-6 })
    {
      const double T0 = Tb + d;
      worst_b = std::fmax( worst_b, std::fabs( ce_p.inverse_integral( ce_p.integral(T0), Tb - 50.0 ) - T0 ) / T0 );
    }
    check( worst_b <= 1e-15, "inverse_integral round trip at the Ce jumps", worst_b, 1e-15 );
  }
  //    c) cost: from the explicit-step guess (as in the Te update kernel), Newton needs a few evaluations even when
  //       the guess is already converged (a sub-ulp step once fell back to ~50 bisections)
  for(double rel : { 0.0, 1e-12, 1e-8, 1e-4, -1e-4, 1e-2, -1e-2 })
  {
    int max_eval = 0;
    for(int i=0;i<=10000;i++)
    {
      const double T0 = 100.0 + 0.99*i;
      const double E1 = ce_p.integral(T0) * ( 1.0 + rel );
      int n = 0;
      ce_p.inverse_integral( E1, T0 + ( E1 - ce_p.integral(T0) ) / ce_p.value(T0), &n );
      max_eval = n > max_eval ? n : max_eval;
    }
    char buf[64]; std::snprintf(buf,sizeof(buf),"inverse_integral evaluations, dE/E=%g",rel);
    check( max_eval <= 4, buf, max_eval, 4 );
  }
  //    d) E <= 0 (energy exhausted) maps to 0 K
  check( ce_p.inverse_integral( 0.0, 300.0 ) == 0.0 && ce_p.inverse_integral( -1.0, 300.0 ) == 0.0, "inverse_integral(E<=0) == 0", 0, 0 );

  std::printf( "%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures==1 ? "" : "s" );
  return g_failures ? 1 : 0;
}
