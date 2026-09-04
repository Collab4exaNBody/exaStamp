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

#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace exaStamp
{
  // Minimal writer for numpy's .npy v1.0 format (magic + version + header dict + raw data) -- only
  // ever need a single flat, uncompressed, homogeneous-dtype array here, not the full npz/zip
  // machinery a general-purpose library (e.g. cnpy, already vendored for the mlip-pace plugin but
  // scoped to a build option we don't want to couple this plugin to) would bring in. Assumes a
  // little-endian host (true for every platform this code targets). Shared within mlip-pod only
  // (write_descriptor_pod.cu, write_descriptor_pod_global.cu) -- SNAP keeps its own independent copy.
  inline void write_npy( const std::string & path, const std::vector<size_t> & shape,
                         const char * descr, const void * data, size_t elem_size )
  {
    size_t nelem = 1; for( size_t s : shape ) nelem *= s;
    std::ostringstream header;
    header << "{'descr': '" << descr << "', 'fortran_order': False, 'shape': (";
    if( shape.size() == 1 ) header << shape[0] << ",";
    else for( size_t i=0; i<shape.size(); i++ ) { header << shape[i]; if( i+1<shape.size() ) header << ", "; }
    header << "), }";
    std::string h = header.str();
    const size_t prefix = 6+2+2; // magic + version + header-length field
    const size_t pad = (64 - (prefix + h.size() + 1) % 64) % 64;
    h.append( pad, ' ' );
    h.push_back( '\n' );

    std::ofstream fout( path, std::ios::binary );
    fout.write( "\x93NUMPY", 6 );
    const char version[2] = { 1, 0 };
    fout.write( version, 2 );
    const uint16_t header_len = static_cast<uint16_t>( h.size() );
    fout.write( reinterpret_cast<const char*>(&header_len), 2 );
    fout.write( h.data(), h.size() );
    fout.write( reinterpret_cast<const char*>(data), nelem*elem_size );
  }
}
