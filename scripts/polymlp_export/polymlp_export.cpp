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

// polymlp_export <polymlp.lammps | polymlp.yaml> <out.yaml>
//
// Exports a polynomial machine-learning potential (A. Seko, pypolymlp / lammps-polymlp-package) for exaStamp's
// polymlp_init: the feature enumeration of the polymlp library (pypolymlp C++ core, BSD-3-Clause, linked by this tool
// only) is written as flat tables, so that exaStamp needs neither the library nor its gtinv data:
//   per centre type t: the order parameters a (n, l, m, type pair; pair features: n, type pair), the features
//   d_k = sum_s c_s prod_q a_{s,q} (gtinv) or d_k = a_k (pair), the polynomial E_t = sum_u w_u prod d_{k(u)}.

#include "polymlp_parse_potential.h"
#include "polymlp_parse_potential_legacy.h"
#include "polymlp_features.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

template<class V> static void flow(std::ostream& o, const V& v)
{
  o << "[";
  for(size_t i=0;i<v.size();i++) { if( i ) o << ", "; o << v[i]; }
  o << "]";
}

int main(int argc, char** argv)
{
  if( argc != 3 ) { std::cerr << "usage: polymlp_export <polymlp file> <out.yaml>" << std::endl; return 1; }
  feature_params fp;
  vector1d coeffs, mass;
  std::vector<std::string> ele;
  try
  {
    if( check_polymlp_legacy( argv[1] ) ) parse_polymlp_legacy( argv[1], fp, coeffs, ele, mass );
    else parse_polymlp( argv[1], fp, coeffs, ele, mass );
  }
  catch( const std::exception& e ) { std::cerr << "polymlp_export: cannot read " << argv[1] << ": " << e.what() << std::endl; return 1; }
  if( fp.pair_type != "gaussian" ) { std::cerr << "polymlp_export: only gaussian radial functions are supported" << std::endl; return 1; }

  Features feat( fp, true );
  Maps& maps = feat.get_maps();
  const bool gtinv = fp.feature_type == "gtinv";

  std::ofstream o( argv[2] );
  o << std::setprecision(17);
  o << "format: exastamp-polymlp\nversion: 1\n";
  o << "source: \"" << argv[1] << "\"\n";
  o << "elements: "; flow( o, ele ); o << "\n";
  o << "masses: "; flow( o, mass ); o << "\n";
  o << "cutoff: " << fp.cutoff << "\n";
  o << "feature_type: " << fp.feature_type << "\n";
  o << "model_type: " << fp.model_type << "\nmax_p: " << fp.maxp << "\nmax_l: " << fp.maxl << "\n";
  o << "gaussians:\n";
  for(const auto& p : fp.params) { o << "  - "; flow( o, p ); o << "\n"; }
  o << "type_pairs:\n";
  for(size_t tp=0;tp<maps.tp_to_types.size();tp++) { o << "  - { types: "; flow( o, maps.tp_to_types[tp] ); o << ", n: "; flow( o, maps.tp_to_n[tp] ); o << " }\n"; }

  size_t nfeat_tot = 0, nterm_tot = 0, npoly_tot = 0;
  o << "types:\n";
  for(int t=0;t<fp.n_type;t++)
  {
    const auto& mt = maps.maps_type[t];
    o << "  - type: " << t << "\n";
    // order parameters: [n, l, m, tp, conj] (gtinv) or [n, tp] (pair)
    std::vector<int> attrs;
    if( gtinv ) for(const auto& a : mt.nlmtp_attrs) { attrs.push_back(a.n); attrs.push_back(a.lm.l); attrs.push_back(a.lm.m); attrs.push_back(a.tp); attrs.push_back(a.lm.conj ? 1 : 0); }
    else for(const auto& a : mt.ntp_attrs) { attrs.push_back(a.n); attrs.push_back(a.tp); }
    o << "    order_parameters: "; flow( o, attrs ); o << "\n";
    // features: CSR over terms, each term: coefficient, ids
    std::vector<long> feat_off{0}, term_off{0};
    std::vector<double> term_coef;
    std::vector<int> term_ids;
    for(const auto& f : mt.features)
    {
      for(const auto& s : f) { term_coef.push_back( s.coeff ); for(int id : s.nlmtp_ids) term_ids.push_back( id ); term_off.push_back( term_ids.size() ); }
      feat_off.push_back( term_coef.size() );
    }
    o << "    feature_offsets: "; flow( o, feat_off ); o << "\n";
    o << "    term_coefficients: "; flow( o, term_coef ); o << "\n";
    o << "    term_offsets: "; flow( o, term_off ); o << "\n";
    o << "    term_ids: "; flow( o, term_ids ); o << "\n";
    // polynomial: coefficient, feature ids. The polymlp readers double the regression coefficients (for their m <= 0
    // half sums): the true ones are written.
    std::vector<double> pc;
    std::vector<long> poff{0};
    std::vector<int> pids;
    for(const auto& u : mt.polynomial) { pc.push_back( 0.5 * coeffs[u.global_id] ); for(int k : u.local_ids) pids.push_back( k ); poff.push_back( pids.size() ); }
    o << "    poly_coefficients: "; flow( o, pc ); o << "\n";
    o << "    poly_offsets: "; flow( o, poff ); o << "\n";
    o << "    poly_features: "; flow( o, pids ); o << "\n";
    nfeat_tot += mt.features.size(); nterm_tot += term_coef.size(); npoly_tot += pc.size();
  }
  std::cout << "polymlp_export: " << ele.size() << " element(s), " << fp.feature_type << ", " << nfeat_tot << " features, " << nterm_tot
            << " feature terms, " << npoly_tot << " polynomial terms (" << coeffs.size() << " coefficients) -> " << argv[2] << std::endl;
  return 0;
}
