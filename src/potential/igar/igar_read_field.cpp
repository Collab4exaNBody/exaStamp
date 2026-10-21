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
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/math/basic_types_stream.h>

#include "igar_field.h"

#include <mpi.h>
#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace exaStamp
{

  using namespace exanb;

  class IgarReadField : public OperatorNode
  {
    ADD_SLOT(MPI_Comm, mpi, INPUT, MPI_COMM_WORLD);
    ADD_SLOT(std::string, file_name, INPUT, REQUIRED, DocString{"Legacy VTK file (STRUCTURED_POINTS, one scalar per point, ASCII or BINARY)."});
    ADD_SLOT(std::string, storage, INPUT, std::string("uint8"), DocString{"Storage of the field values: uint8 (1 byte per value, for 8-bit images: every value of the file must be an integer 0-255) or double."});
    ADD_SLOT(double, value_scale, INPUT, 1.0, DocString{"The field is the file values times value_scale (e.g. 1/255 to map the levels 0-255 of an 8-bit image onto [0,1])."});
    ADD_SLOT(std::string, grid_ordering, INPUT, std::string("F_ORDER"), DocString{"Ordering of the values in the file: F_ORDER (x fastest) or C_ORDER (z fastest)."});
    ADD_SLOT(IgarField, igar_field, INPUT_OUTPUT, DocString{"IGAR field, replicated on every MPI process."});

    // values are passed to store(file_index, value) in file order
    template <class StoreT> static inline void parse_ascii(std::ifstream &file, size_t count, StoreT &store, const std::string &fname)
    {
      static constexpr size_t chunk = size_t(1) << 26;
      std::vector<char> buf(chunk + 64);
      size_t keep = 0, n = 0;
      bool eof = false;
      while (n < count && !eof)
      {
        file.read(buf.data() + keep, chunk - keep);
        const size_t len = keep + size_t(file.gcount());
        eof = (file.gcount() == 0) || file.eof();
        const char *p = buf.data();
        const char *end = buf.data() + len;
        while (n < count)
        {
          while (p < end && (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t')) ++p;
          // a number may be cut at the chunk end: keep it for the next chunk
          const char *q = p;
          while (q < end && !(*q == ' ' || *q == '\n' || *q == '\r' || *q == '\t')) ++q;
          if (q == end && !eof) break;
          if (p == end) break;
          double v = 0.0;
          const auto res = std::from_chars(p, q, v);
          if (res.ec != std::errc() || res.ptr != q) fatal_error() << "igar_read_field: " << fname << ": bad value '" << std::string(p, q) << "' at index " << n << std::endl;
          store(n++, v);
          p = q;
        }
        keep = size_t(end - p);
        std::memmove(buf.data(), p, keep);
      }
      if (n != count) fatal_error() << "igar_read_field: " << fname << ": " << n << " values read, " << count << " expected" << std::endl;
    }

    template <class T, class StoreT> static inline void parse_binary_typed(std::ifstream &file, size_t count, StoreT &store, const std::string &fname)
    {
      static constexpr size_t chunk_values = (size_t(1) << 26) / sizeof(T);
      std::vector<unsigned char> buf(chunk_values * sizeof(T));
      size_t n = 0;
      while (n < count)
      {
        const size_t nv = std::min(chunk_values, count - n);
        file.read(reinterpret_cast<char *>(buf.data()), nv * sizeof(T));
        if (size_t(file.gcount()) != nv * sizeof(T)) fatal_error() << "igar_read_field: " << fname << ": unexpected end of file after " << n << " values" << std::endl;
        for (size_t i = 0; i < nv; i++)
        {
          // legacy VTK binary data is big-endian
          unsigned char be[sizeof(T)];
          for (size_t b = 0; b < sizeof(T); b++) be[b] = buf[i * sizeof(T) + sizeof(T) - 1 - b];
          T x;
          std::memcpy(&x, be, sizeof(T));
          store(n + i, static_cast<double>(x));
        }
        n += nv;
      }
    }

    template <class StoreT> static inline void parse_binary(std::ifstream &file, const std::string &type, size_t count, StoreT &store, const std::string &fname)
    {
      if (type == "unsigned_char" || type == "vtktypeuint8") parse_binary_typed<uint8_t>(file, count, store, fname);
      else if (type == "char" || type == "vtktypeint8") parse_binary_typed<int8_t>(file, count, store, fname);
      else if (type == "unsigned_short" || type == "vtktypeuint16") parse_binary_typed<uint16_t>(file, count, store, fname);
      else if (type == "short" || type == "vtktypeint16") parse_binary_typed<int16_t>(file, count, store, fname);
      else if (type == "unsigned_int" || type == "vtktypeuint32") parse_binary_typed<uint32_t>(file, count, store, fname);
      else if (type == "int" || type == "vtktypeint32") parse_binary_typed<int32_t>(file, count, store, fname);
      else if (type == "vtktypeuint64" || type == "unsigned_long") parse_binary_typed<uint64_t>(file, count, store, fname);
      else if (type == "vtktypeint64" || type == "long") parse_binary_typed<int64_t>(file, count, store, fname);
      else if (type == "float" || type == "vtktypefloat32") parse_binary_typed<float>(file, count, store, fname);
      else if (type == "double" || type == "vtktypefloat64") parse_binary_typed<double>(file, count, store, fname);
      else fatal_error() << "igar_read_field: " << fname << ": unsupported binary data type '" << type << "'" << std::endl;
    }

    // rank 0 only: fills the field from the file
    inline void read_file(IgarField &field, bool as_uint8, bool c_order)
    {
      const std::string &fname = *file_name;
      std::ifstream file(fname, std::ios::binary);
      if (!file.is_open()) fatal_error() << "igar_read_field: could not open " << fname << std::endl;

      bool binary = false, lookup_table = false, color_scalars = false;
      IJK dims = {0, 0, 0};
      size_t count = 0;
      std::string type;
      std::string line;
      while (std::getline(file, line))
      {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::istringstream ls(line);
        std::string keyword;
        ls >> keyword;
        if (keyword == "BINARY") binary = true;
        else if (keyword == "DATASET")
        {
          std::string ds;
          ls >> ds;
          if (ds != "STRUCTURED_POINTS") fatal_error() << "igar_read_field: " << fname << ": DATASET " << ds << " not supported (STRUCTURED_POINTS only)" << std::endl;
        }
        else if (keyword == "DIMENSIONS") ls >> dims.i >> dims.j >> dims.k;
        else if (keyword == "POINT_DATA") ls >> count;
        else if (keyword == "SCALARS")
        {
          std::string name;
          int ncomp = 1;
          ls >> name >> type;
          if (ls >> ncomp && ncomp != 1) fatal_error() << "igar_read_field: " << fname << ": scalar field expected, " << ncomp << " components found" << std::endl;
        }
        else if (keyword == "LOOKUP_TABLE")
        {
          lookup_table = true;
          break;
        }
        else if (keyword == "COLOR_SCALARS")
        {
          // 8-bit data written by VTK for unsigned char arrays: bytes in BINARY files, byte/255 in ASCII files, no LOOKUP_TABLE
          std::string name;
          int ncomp = 0;
          ls >> name >> ncomp;
          if (ncomp != 1) fatal_error() << "igar_read_field: " << fname << ": scalar field expected, COLOR_SCALARS with " << ncomp << " values found" << std::endl;
          type = "unsigned_char";
          color_scalars = true;
          lookup_table = true;
          break;
        }
      }
      if (!lookup_table) fatal_error() << "igar_read_field: " << fname << ": LOOKUP_TABLE not found" << std::endl;
      if (dims.i <= 0 || dims.j <= 0 || dims.k <= 0) fatal_error() << "igar_read_field: " << fname << ": bad DIMENSIONS " << dims << std::endl;
      field.m_dims = dims;
      if (count != field.number_of_values()) fatal_error() << "igar_read_field: " << fname << ": POINT_DATA " << count << " != " << dims.i << "x" << dims.j << "x" << dims.k << std::endl;

      if (as_uint8) field.m_u8.resize(count);
      else field.m_f64.resize(count);
      uint8_t *const u8 = field.m_u8.data();
      double *const f64 = field.m_f64.data();
      const size_t nx = dims.i, ny = dims.j, nz = dims.k;

      auto store = [&](size_t n, double v)
      {
        size_t dest = n;
        if (c_order) dest = n / (ny * nz) + nx * ((n / nz) % ny + ny * (n % nz));
        if (as_uint8)
        {
          const double b = std::nearbyint(v);
          if (!(b >= 0.0 && b <= 255.0 && std::fabs(b - v) <= 1.0e-6))
          {
            fatal_error() << "igar_read_field: " << fname << ": value " << v << " at index " << n << " is not an integer 0-255 (use storage: double for a non 8-bit field)" << std::endl;
          }
          u8[dest] = static_cast<uint8_t>(b);
        }
        else
        {
          f64[dest] = v;
        }
      };

      if (binary) parse_binary(file, type, count, store, fname);
      else if (color_scalars)
      {
        auto store_color = [&](size_t n, double v) { store(n, std::nearbyint(v * 255.0)); };
        parse_ascii(file, count, store_color, fname);
      }
      else parse_ascii(file, count, store, fname);
    }

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Reads the IGAR scalar field used by igar_force_interp from a legacy VTK file
(DATASET STRUCTURED_POINTS, one scalar per point, ASCII or BINARY; SCALARS of
any numeric type, or COLOR_SCALARS with one value per point, read as bytes).

The field covers the whole simulation domain: its DIMENSIONS nx ny nz split
the domain bounds into nx*ny*nz voxels, each value sitting at a voxel center
(SPACING and ORIGIN of the file are ignored). Process 0 reads the file and
broadcasts it: every MPI process holds the whole field, which is never
exchanged nor migrated afterwards.

storage selects how values are kept in memory:
  uint8 (default)  1 byte per value. For 8-bit images: every value of the file
                   must be an integer 0-255, otherwise the run stops (values
                   are never rounded). Energies and forces are the same as
                   with double storage. A 1024^3 field takes 1 GiB per process.
  double           8 bytes per value, any real values.

The field is the file values times value_scale (default 1). An image stored
normalized to [0,1] can be turned into bytes with vtk_ascii_to_binary.py
--scale 255, then read with storage: uint8, value_scale: 1/255.

Example:

  setup_system:
    - igar_read_field:
        file_name: "image_1024.vtk"
        storage: uint8
)EOF";
    }

    inline void execute() override final
    {
      MPI_Comm comm = *mpi;
      int rank = 0;
      MPI_Comm_rank(comm, &rank);

      bool as_uint8 = true;
      if (*storage == "uint8") as_uint8 = true;
      else if (*storage == "double") as_uint8 = false;
      else fatal_error() << "igar_read_field: unknown storage '" << *storage << "' (uint8 or double)" << std::endl;
      bool c_order = false;
      if (*grid_ordering == "C_ORDER") c_order = true;
      else if (*grid_ordering != "F_ORDER") fatal_error() << "igar_read_field: unknown grid_ordering '" << *grid_ordering << "' (F_ORDER or C_ORDER)" << std::endl;
      const double scale = *value_scale;
      if (!(scale > 0.0)) fatal_error() << "igar_read_field: value_scale must be > 0" << std::endl;

      IgarField &field = *igar_field;
      field = IgarField{};
      field.m_value_scale = scale;

      const auto t0 = std::chrono::steady_clock::now();
      if (rank == 0) read_file(field, as_uint8, c_order);

      long long dims[3] = {field.m_dims.i, field.m_dims.j, field.m_dims.k};
      MPI_Bcast(dims, 3, MPI_LONG_LONG, 0, comm);
      field.m_dims = IJK{dims[0], dims[1], dims[2]};
      const size_t count = field.number_of_values();
      unsigned char *data = nullptr;
      size_t bytes = 0;
      if (as_uint8)
      {
        field.m_u8.resize(count);
        data = field.m_u8.data();
        bytes = count;
      }
      else
      {
        field.m_f64.resize(count);
        data = reinterpret_cast<unsigned char *>(field.m_f64.data());
        bytes = count * sizeof(double);
      }
      static constexpr size_t bcast_chunk = size_t(1) << 28;
      for (size_t off = 0; off < bytes; off += bcast_chunk)
      {
        MPI_Bcast(data + off, int(std::min(bcast_chunk, bytes - off)), MPI_BYTE, 0, comm);
      }
      const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

      lout << "IGAR field       = " << *file_name << std::endl
           << "IGAR dimensions  = " << field.m_dims.i << " x " << field.m_dims.j << " x " << field.m_dims.k << std::endl
           << "IGAR storage     = " << *storage << " (" << (bytes >> 20) << " MiB per process), value_scale = " << scale << std::endl
           << "IGAR read time   = " << elapsed << " s" << std::endl;
    }
  };

  ONIKA_AUTORUN_INIT(igar_read_field) { OperatorNodeFactory::instance()->register_factory("igar_read_field", make_simple_operator<IgarReadField>); }

} // namespace exaStamp
