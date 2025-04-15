#pragma once

#include <array>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <utility>

#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/string_utils.h>
#include <onika/file_utils.h>
#include <onika/log.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>

#include <exaStamp/io/tokenizer_utils.h>

#define IOEXTRA_FILE_BUFFER_SIZE 8192
#define IOEXTRA_NULL_CHAR '\0'
#define IOEXTRA_LZMA_STREAM_BUF_SIZE 8192
#define IOEXTRA_LZMA_INTERNAL_BUF_SIZE 8192
#define IOEXTRA_BZ2_STREAM_BUF_SIZE 8192
#define IOEXTRA_BZ2_INTERNAL_BUF_SIZE 8192

// helper macro for debugging

#define lerr()                                                                                                         \
  onika::lerr << onika::format_string("[%s:%s:%s] ", std::filesystem::path(__FILE__).filename().c_str(), __func__,     \
                                      std::to_string(__LINE__).c_str())

#define lwarn() onika::lout << "[WARNING] "

namespace exaStamp {

using namespace exanb;


// using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_id, field::_type>;
using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_vx, field::_vy, field::_vy, field::_vz, field::_id, field::_type>;
using ParticleData = std::vector<ParticleTupleIO>;

namespace ioextra {

using namespace tokens; // tokenizer utilities

/* ------------------------------------------------------------------------- */

template <size_t N> struct SpeciesTmpl : TokenSetTmpl<N> {
  inline bool get_type_from_sv(std::string_view type, size_t& index) {
    return match_token_any(type, this->tokens, index);
  }
};

using Species = SpeciesTmpl<16>;

// ExaStamp context when parsing a file.
struct Context {

  ParticleData& particle_data;

  Species species;

  size_t n_particles = 0;
  size_t n_species = 0; // number of species detected in the file

  Mat3d H = make_identity_matrix(); // lattice vector
  Vec3d origin{0., 0., 0.};         // box origin

  Mat3d Hinv{};
  Mat3d D{};
  Mat3d D2{};
  Mat3d Dinv{};
  
  Mat3d Xform = Mat3d{};            // initial xform (likely == I)
  Mat3d invXform = Mat3d{};         // intial xform^-1
  Vec3d boxlen{};                   // box len in x, y, z

  bool uniform_scale;               //

  // for timing purpose
  size_t elapsed_time_ms = 0;
  size_t elapsed_time_ys = 0;
  size_t elapsed_time_ns = 0;

  inline void init_domain(Domain& domain, ReadBoundsSelectionMode bounds_mode) {

    // assume H follow column vector convention
    Vec3d a{H.m11, H.m21, H.m31};
    Vec3d b{H.m12, H.m22, H.m32};
    Vec3d c{H.m13, H.m23, H.m33};

    boxlen.x = norm(a);
    boxlen.y = norm(b);
    boxlen.z = norm(c);

    lout << "a = " << a << onika::format_string(" , norm = %10.5f ", boxlen.x) << std::endl;
    lout << "b = " << b << onika::format_string(" , norm = %10.5f ", boxlen.y) << std::endl;
    lout << "c = " << c << onika::format_string(" , norm = %10.5f ", boxlen.z) << std::endl;

    D = diag_matrix(boxlen);
    Dinv = inverse(D);
    Mat3d F1 = H * Dinv;

    double cellsize = domain.cell_size();

    // check for grid dims. If grid dims is not provided by the use, 
    // it will be deducde from box lengths and cellsize.
    IJK grid_dim = domain.grid_dimension();
    size_t nx = static_cast<ssize_t>(boxlen.x / cellsize);
    size_t ny = static_cast<ssize_t>(boxlen.y / cellsize);
    size_t nz = static_cast<ssize_t>(boxlen.z / cellsize);
    nx = (grid_dim.i > 0 && cmp::ne(grid_dim.i, nx)) ? grid_dim.i : nx;
    ny = (grid_dim.j > 0 && cmp::ne(grid_dim.j, ny)) ? grid_dim.j : ny;
    nz = (grid_dim.k > 0 && cmp::ne(grid_dim.k, nz)) ? grid_dim.k : nz;

    D2 = diag_matrix(Vec3d{nx * cellsize, ny * cellsize, nz * cellsize});
    Mat3d D2inv = inverse(D2);
    Mat3d F2 = D * D2inv;

    Xform = F1 * F2;
    invXform = inverse(Xform);

    AABB domain_bounds = {{0., 0., 0.}, {D2.m11, D2.m22, D2.m33}};
    compute_domain_bounds(domain, bounds_mode, 0.0, domain_bounds, domain_bounds, true);
    uniform_scale = is_uniform_scale(Xform);

    if (uniform_scale) {
      domain.set_xform(make_identity_matrix());
      domain.set_cell_size(domain.cell_size() * Xform.m11);
      domain.set_bounds({domain.origin() * Xform.m11, domain.extent() * Xform.m11});
    } else {
      domain.set_xform(Xform);
    }
  }
};

using IOContext = Context;

inline void wrap_to_domain(Vec3d& r, const IOContext& ctx) {
    r.x = wrap(r.x, ctx.D2.m11);
    r.y = wrap(r.y, ctx.D2.m22);
    r.z = wrap(r.z, ctx.D2.m33);
}

/* ------------------------------------------------------------------------- */

// mode to open a file
enum FileMode : char {
  READ   = 'r',
  WRITE  = 'w',
  APPEND = 'a',
};

// compression mode
enum FileCompression {
  NONE,
  GZIP,
  BZIP2,
  XZ,
};

inline const char* convert_mode_to_char(FileMode mode) {
  switch (mode) {
  case FileMode::READ:
    return "rb";
    break;
  case FileMode::APPEND:
    return "a+b";
    break;
  case FileMode::WRITE:
    return "wb";
    break;
  default:
    lerr << "Invalid convertion from Filemode to char *" << std::endl;
    std::abort();
    break;
  }
}

inline FileCompression convert_char_to_compression(const std::string& str) {
  if (str == "gz")
    return FileCompression::GZIP;
  if (str == "bz2")
    return FileCompression::BZIP2;
  if (str == "xz")
    return FileCompression::XZ;
  return FileCompression::NONE;
}

inline std::string convert_compression_to_char(FileCompression& compression) {
  if (compression == FileCompression::GZIP)
    return "gz";
  if (compression == FileCompression::BZIP2)
    return "bz2";
  if (compression == FileCompression::XZ)
    return "xz";
  return "";
}

/* ------------------------------------------------------------------------- */

// Base class that handler basic file operation
// Specialization allow abstraction on file types.
class TextFileHandler {
public:
  TextFileHandler(const std::string& filepath);

  virtual ~TextFileHandler() = default;

  virtual void clear() noexcept = 0;
  virtual void seek(uint64_t position) = 0;
  virtual size_t read(char* buffer, size_t size) = 0;

protected:
  const std::string& path() { return m_path; };
  std::string m_path;
};

class ASCIIFileHandler final : public TextFileHandler {
public:
  ASCIIFileHandler(const std::string& path, FileMode mode);
  ~ASCIIFileHandler() override;

  void clear() noexcept override;
  void seek(uint64_t cursor) override;

  size_t read(char* data, size_t size) override;

private:
  std::FILE* m_fptr;
};


#ifdef USE_ZLIB
#define ZLIB_CONST
#include <zconf.h>
#include <zlib.h>

typedef struct gzFile_s *gzFile;

class GzipFileHandler final : public TextFileHandler {
public:
  GzipFileHandler(const std::string& path, FileMode mode);
  ~GzipFileHandler() override;

  void clear() noexcept override;
  void seek(uint64_t cursor) override;

  size_t read(char* buffer, size_t size) override;

  const char* gz_error() const;
  static unsigned safe_cast(size_t);
  

private:
  gzFile m_fptr = nullptr;
};

#endif

#ifdef USE_LZMA
#include <lzma.h>

class XzFileHandler final : public TextFileHandler {
public:
  XzFileHandler(const std::string& path, FileMode mode);
  ~XzFileHandler() override;

  void clear() noexcept override;
  void seek(uint64_t cursor) override;

  size_t read(char* buffer, size_t size) override;

  size_t safe_cast(uint64_t value);
  void check_lzma_ret(lzma_ret ret);
  void start_lzma_decoder_stream(lzma_stream* stream);

private:
  std::FILE* m_fptr = nullptr;
  lzma_stream m_xz_stream = LZMA_STREAM_INIT;
  
  std::vector<uint8_t> m_xz_buffer = {0};

  // static constexpr size_t m_xz_buffer_end = IOEXTRA_LZMA_INTERNAL_BUF_SIZE;
  // uint8_t* m_xz_buffer[m_xz_buffer_end];
};

#endif

#ifdef USE_BZIP2
#include <bzlib.h>

class Bzip2FileHandler final : public TextFileHandler {
public:
  Bzip2FileHandler(const std::string& path, FileMode mode);
  ~Bzip2FileHandler() override;

  void clear() noexcept override;
  void seek(uint64_t cursor) override;

  size_t read(char* buffer, size_t size) override;

  unsigned safe_cast(uint64_t size);
  void check_bz2_retcode(int code);

private:
  std::FILE* m_fptr;
  std::function<int(bz_stream*)> m_end_bz2_stream;
  bz_stream m_bz2_stream;
  std::vector<char> m_bz2_buffer;
};

#endif

/* ------------------------------------------------------------------------- */

class File {
protected:
  std::string m_path;
  FileMode m_mode;
  FileCompression m_compression;

  File(std::string path, FileMode mode, FileCompression compression)
      : m_path(std::move(path)), m_mode(mode), m_compression(compression) {}

public:
  inline const std::string& path() const { return m_path; };
  inline const FileMode mode() const { return m_mode; };
  inline const FileCompression compression() const { return m_compression; };
};

class TextFile : public File {
private:
  std::vector<char> m_buf;
  const char* m_lstart; // line start
  const char* m_bend;   // buf end

  uint64_t m_cursor = 0;
  bool m_eof = false;
  bool m_hdlr_eof = false;

  std::unique_ptr<TextFileHandler> m_handler;

  inline bool is_buffer_init() const { return m_buf[0] != IOEXTRA_NULL_CHAR; }
  void fill_buffer(size_t pos);

public:
  TextFile(std::string filepath, FileMode mode, FileCompression compression);

  inline bool eof() const { return m_eof; };

  uint64_t tell() const;
  void seek(uint64_t pos);
  void clear();
  void reset();

  // get the next line in the buffer.
  // If the line is not complete, load the next section of the buffer
  std::string_view get_line();
};

/* ------------------------------------------------------------------------- */

struct Metadata {
  const std::string name = "";
  const std::vector<std::string> aliases = {};
  const std::string description = "";
};

template<typename T>
const Metadata get_parser_metadata() { 
  return { };
}

// Base class to all format parsers
class Parser {
public:
  Parser() = default;
  virtual ~Parser() = default;

  Parser(const Parser&) = delete;
  Parser& operator=(const Parser&) = delete;
  Parser(Parser&&) = delete;
  Parser& operator=(Parser&&) = delete;

  // read the first step in the file
  virtual inline bool operator()(Context& ctx) { return false; };
  // read the step at the given index
  virtual inline bool operator()(Context& ctx, size_t index) { return false; };
  // returns the number of known step in the file
  virtual const size_t size() = 0;
};

// Specialization for format that are just text files
class TextParser : public Parser {

protected:
  bool m_eof = false;
  std::vector<size_t> loc{};
  TextFile m_file;

  std::array<std::string_view, 10> line_buffer{};
  inline std::string_view& current_line(void) { return line_buffer[0]; };

  std::string_view m_current_line;
  TokenSet m_tokens{};

public:
  TextParser(std::string filepath, FileMode mode, FileCompression compression)
      : m_file(std::move(filepath), mode, compression) {}

  virtual ~TextParser() override = default;

  virtual int64_t next() = 0;
  virtual bool parse(Context& ctx) = 0;

  bool operator()(Context& ctx) override;
  bool operator()(Context& ctx, size_t index) override;

  void scan();

  inline const size_t size() override { return loc.size(); };
  inline const std::string_view& current_line(void) const { return line_buffer[0]; }
};

/* ------------------------------------------------------------------------- */

using ParserCreator = std::function<std::unique_ptr<Parser>(std::string, FileMode, FileCompression)>;

struct RegisteredFormat {
  const Metadata metadata;
  const ParserCreator creator;
  bool is_valid = true;
};

class ParserFactory {
public:

  ParserFactory();

  std::unordered_map<std::string, size_t> alias_to_index;
  std::unordered_map<size_t, RegisteredFormat> index_to_parser;

  template <typename T> void register_format() {

    size_t index = index_to_parser.size();
    Metadata metadata = get_parser_metadata<T>();

    index_to_parser.insert(
        {index, {.metadata = metadata, .creator = [](std::string filepath, FileMode mode, FileCompression compression) {
                   return std::make_unique<T>(filepath, mode, compression);
                 }}});

    // map extension to the index
    for (const std::string& alias : metadata.aliases) {
      std::unordered_map<std::string, size_t>::const_iterator it = alias_to_index.find(alias);
      if (it != alias_to_index.end()) {
        continue; // the alias is already used
      }
      alias_to_index.insert({alias, index});
    }
  }

  RegisteredFormat from_alias(const std::string& alias) const {
    std::unordered_map<std::string, size_t>::const_iterator it = alias_to_index.find(alias);
    if (it == alias_to_index.end()) {
      return {.is_valid = false};
    }
    return index_to_parser.find(it->second)->second;
  }

  inline const bool has_alias(const std::string& alias) const {
    return alias_to_index.find(alias) != alias_to_index.end();
  }

  std::string infos() const;
};

const ParserFactory& parser_factory();

/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Data format

namespace lmpdata {

enum class AtomStyle { ATOMIC, MOLECULAR, FULL, NONE };

enum class Field { NATOM, TYPAT, XAXIS, YAXIS, ZAXIS, TILT, SECTION, END };

enum class Section {
  HEADER,
  MASSES,
  ATOMS,
  VELOCITIES,
  IGNORED,
};

constexpr std::array<Section, 3> sections_list { Section::MASSES, Section::ATOMS, Section::VELOCITIES };

constexpr TokenNeedles<1> tok_atoms = {"atoms"};
constexpr TokenNeedles<2> tok_atom_types = {"atom", "types"};
constexpr TokenNeedles<2> tok_xaxis = {"xlo", "xhi"};
constexpr TokenNeedles<2> tok_yaxis = {"ylo", "yhi"};
constexpr TokenNeedles<2> tok_zaxis = {"zlo", "zhi"};
constexpr TokenNeedles<3> tok_tilt = {"xy", "xz", "yz"};
constexpr TokenNeedles<3> tok_sections = {"Masses", "Atoms", "Velocities"};
constexpr TokenNeedles<3> tok_atom_style = {"atomic", "full", "molecule"};


class LAMMPSDataParser : public TextParser {
public:
  LAMMPSDataParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression), m_atom_style(AtomStyle::NONE),
        m_current_section(Section::HEADER) {}

  int64_t next() override;
  bool parse(Context& proxy) override;

private:

  AtomStyle m_atom_style;
  Section m_current_section;
  
  // keep track of the current section line
  inline std::string_view& section_line(void) { return line_buffer[1]; };
  inline const std::string_view& section_line(void) const { return line_buffer[1]; };

  bool parse_header(Context& ctx);
  bool parse_atoms(Context& ctx);
  bool parse_masses();
  bool parse_velocities(Context& ctx);

  void forward_to_next_section();
};

using AtomLineReader = bool(*)(const TokenSet&, Context&, size_t);
using FieldParser = bool(*)(const TokenSet&, Context&);

bool read_atom_line_atomic(const TokenSet&, Context&, size_t);

inline const std::unordered_map<AtomStyle, AtomLineReader>& atom_line_readers() {
  static std::unordered_map<AtomStyle, AtomLineReader> map = {
    { AtomStyle::ATOMIC , &read_atom_line_atomic }
  };
  return map;
}

bool read_field_natom(const TokenSet&, Context&);
bool read_field_atom_types(const TokenSet&, Context&);
bool read_field_xaxis(const TokenSet&, Context&);
bool read_field_yaxis(const TokenSet&, Context&);
bool read_field_zaxis(const TokenSet&, Context&);
bool read_field_tilt(const TokenSet&, Context&);
bool read_field_section(const TokenSet&, Context&);

inline const std::array<std::pair<Field, FieldParser>, 7>& field_readers() {
  static constexpr std::array<std::pair<Field, FieldParser>, 7> readers = {{
      {Field::NATOM, &read_field_natom},
      {Field::TYPAT, &read_field_atom_types},
      {Field::XAXIS, &read_field_xaxis},
      {Field::YAXIS, &read_field_zaxis},
      {Field::ZAXIS, &read_field_yaxis},
      {Field::TILT, &read_field_tilt},
      {Field::SECTION, &read_field_section},
  }};
  return readers;
}

inline Section get_section(std::string_view str) {
  size_t index;
  if (match_token_any(str, tok_sections, index))
    return sections_list[index];
  else
    return Section::IGNORED;
}

} // namespace lmpdata

template<> inline const Metadata get_parser_metadata<lmpdata::LAMMPSDataParser>() {
  return {
    .name        = "LAMMPS Data",
    .aliases     = { "lmp", "lmp-data" },
    .description = "LAMMPS Data format"
  };
}

/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Dump format
// namespace lmpdump {
// using namespace re;

// static const std::string str_item_natom = std::string("^(?!#)ITEM:\\sNUMBER\\sOF\\sATOMS");
// static const std::string str_item_bbox = std::string("^(?!#)ITEM:\\sBOX\\sBOUNDS");
// static const std::string str_item_atoms = std::string("^(?!#)ITEM:\\sATOMS");
// static const std::string str_item_tilt = std::string("xy\\sxz\\syz");
// static const std::string str_item_fields = std::string("([\\w\\[\\]]+)");

// const std::regex re_item_natom = std::regex(str_item_natom);
// const std::regex re_item_bbox  = std::regex(str_item_bbox);
// const std::regex re_item_tilt  = std::regex(str_item_tilt);
// const std::regex re_item_atoms = std::regex(str_item_atoms);
// const std::regex re_item_fields = std::regex(str_item_fields);

// enum class X : uint8_t {
//   NONE = 0,
//   X    = 1 << 0,   // wrapped
//   XS   = 1 << 1,  // scaled ,
//   XU   = 1 << 2,  // unsrapped,
//   XSU  = 1 << 3, // scaled unwrapped
// };

// inline constexpr X operator|(X lhs, X rhs) {
//   return static_cast<X>(static_cast<uint8_t>(lhs) | static_cast<uint8_t>(rhs));
// }

// inline constexpr X operator&(X lhs, X rhs) {
//   return static_cast<X>(static_cast<uint8_t>(lhs) & static_cast<uint8_t>(rhs));
// }

// inline constexpr X operator~(X f) { return static_cast<X>(~static_cast<uint8_t>(f)); }

// inline constexpr X& operator|=(X& lhs, X rhs) {
//   lhs = lhs | rhs;
//   return lhs;
// }

// inline constexpr X& operator&=(X& lhs, X rhs) {
//   lhs = lhs & rhs;
//   return lhs;
// }

// enum class Field {
//   ID,      // atom id
//   TYPE,    // atom type
//   ELEMENT, // element
//   RX,      // x wrapped
//   RY,      // y wrapped
//   RZ,      // z wrapped
//   RXU,     // x unwrapped
//   RYU,     // y unwrapped
//   RZU,     // z unwrapped
//   RXS,     // x scaled position
//   RYS,     // y scaled position
//   RZS,     // z scaled position
//   RXSU,    // x scaled unwrapped position
//   RYSU,    // y scaled unwrapped position
//   RZSU,    // z scaled unwrapped position
//   VX,      // x velocity
//   VY,      // y velocity
//   VZ,      // z velocity
//   IGNORED, // ignored field
// };


// constexpr std::array<std::pair<std::array<Field, 3>, X>, 4> x_priority_groups {{
//   {{ Field::RX  , Field::RY  , Field::RZ  }, X::X  },   // 1st prioprity
//   {{ Field::RXU , Field::RYU , Field::RZU }, X::XU },   // 2nd priority
//   {{ Field::RXS , Field::RYS , Field::RZS }, X::XS },   // 3rd priority
//   {{ Field::RXSU, Field::RYSU, Field::RZSU}, X::XSU},   // 4th priority
// }};

// using FieldMap = std::unordered_map<std::string, Field>;

// inline const FieldMap& field_map() {
//   static const FieldMap map =  {
//     { "id"     , Field::ID      },
//     { "type"   , Field::TYPE    },
//     { "element", Field::ELEMENT },
//     { "x"      , Field::RX      },
//     { "y"      , Field::RY      },
//     { "z"      , Field::RZ      },
//     { "xu"     , Field::RXU     },
//     { "yu"     , Field::RYU     },
//     { "zu"     , Field::RZU     },
//     { "xs"     , Field::RXS     },
//     { "ys"     , Field::RYS     },
//     { "zs"     , Field::RZS     },
//     { "xsu"    , Field::RXSU    },
//     { "ysu"    , Field::RYSU    },
//     { "zsu"    , Field::RZSU    },
//     { "vx"     , Field::VX      },
//     { "vx"     , Field::VX      },
//     { "vx"     , Field::VX      },
//     { "nil"    , Field::IGNORED }
//   };
//   return map;
// };

// struct AtomFields {
//   size_t size = 0;
//   X x = X::NONE;
//   std::array<Field, 3> xfields = {Field::RX, Field::RY, Field::RZ};
//   std::unordered_map<Field, size_t> indices;

//   bool has_type = false;
//   bool has_id = false;
//   bool scaled = false;
//   bool unwrapped = false;
//   bool triclinic = false;

//   inline const Field& xfield() { return xfields[0]; };
//   inline const Field& yfield() { return xfields[1]; };
//   inline const Field& zfield() { return xfields[2]; };
// };

// class LAMMPSDumpParser : public TextParser {
// public:
//   LAMMPSDumpParser(std::string filepath, FileMode mode, FileCompression compression)
//       : TextParser(filepath, mode, compression) {}

//   int64_t next() override;
//   bool parse(Context& ctx) override;

// private:
//   void xyz_fields(double&, double&, double&, const AtomFields&, const Context&);
//   bool process_fields(AtomFields&, const std::string_view&);

//   inline std::string build_line_regex_pattern(size_t size) {
//     std::string pattern = "^\\s*";
//     for (size_t i = 0; i < size; ++i) {
//       pattern += str_any;
//       if (i < size - 1)
//         pattern += str_space;
//     }
//     return pattern;
//   }
// };

// } // namespace lmpdump

// template<> inline const Metadata get_parser_metadata<lmpdump::LAMMPSDumpParser>() {
//   return {
//     .name        = "LAMMPS Dump",
//     .aliases     = { "dump", "lmp-dump", "lammps-dump" },
//     .description = "LAMMPS Dump format"
//   };
// }

/* ------------------------------------------------------------------------- */
// Extended XYZ Format parser

namespace xyz {

enum class FieldType : uint8_t {
  N = 0,
  R = 1 << 0,
  S = 1 << 1,
  I = 1 << 2
};

inline constexpr FieldType operator|(FieldType lhs, FieldType rhs) {
  return static_cast<FieldType>(static_cast<uint8_t>(lhs) | static_cast<uint8_t>(rhs));
}

inline constexpr FieldType operator&(FieldType lhs, FieldType rhs) {
  return static_cast<FieldType>(static_cast<uint8_t>(lhs) & static_cast<uint8_t>(rhs));
}

inline constexpr FieldType operator~(FieldType f) { return static_cast<FieldType>(~static_cast<uint8_t>(f)); }

inline constexpr FieldType& operator|=(FieldType& lhs, FieldType rhs) {
  lhs = lhs | rhs;
  return lhs;
}

inline constexpr FieldType& operator&=(FieldType& lhs, FieldType rhs) {
  lhs = lhs & rhs;
  return lhs;
}

enum class Field { SPECIES, POSITIONS, VELOCITIES, ID, TYPE, NIL };

struct AtomicProperty {
  Field field = Field::NIL;
  FieldType type = FieldType::N;
  size_t size = 0;
  size_t index = 0;
};

struct Properties {
  size_t index = 0;
  std::array<AtomicProperty, 4> data{};

  constexpr inline const AtomicProperty& species() const { return *(data.data() + 0); }
  constexpr inline const AtomicProperty& positions() const { return *(data.data() + 1); }
  constexpr inline const AtomicProperty& velo() const { return *(data.data() + 2); }
  constexpr inline const AtomicProperty& id() const { return *(data.data() + 3); }

  constexpr inline const bool has_species() const { return species().field == Field::SPECIES; }
  constexpr inline const bool has_positions() const { return positions().field == Field::POSITIONS; }
  constexpr inline const bool has_velo() const { return velo().field == Field::VELOCITIES; }
  constexpr inline const bool has_id() const { return positions().field == Field::ID; }
};

// order matter - same order as in Properties
constexpr TokenNeedles<4> tok_props{"species", "pos", "velo", "id"};
constexpr std::array<Field, 4> props_list{Field::SPECIES, Field::POSITIONS, Field::VELOCITIES, Field::ID};

constexpr TokenNeedles<3> tok_field_type = {"R", "S", "I"};
constexpr std::array<FieldType, 3> type_list = {FieldType::R, FieldType::S, FieldType::I};

template <size_t N1, size_t N2>
bool read_exyz_comment_line(std::string_view line, TokenSetTmpl<N1>& tok_lattice, TokenSetTmpl<N2>& tok_props) {
  constexpr std::string_view lattice_key = "Lattice=\"";
  constexpr size_t lattice_key_len = 9;
  constexpr std::string_view prop_key = "Properties=";
  constexpr size_t prop_key_len = 11;

  const char* begin = line.data();
  const char* end = line.data() + line.size();

  if (!starts_with(line, lattice_key))
    return false;

  const char* val_start = begin + lattice_key_len;
  const char* val_end = val_start;

  while (val_end < end && *val_end != '"')
    ++val_end;
  if (val_end >= end)
    return false;

  std::string_view lattice_view = std::string_view(val_start, val_end);
  tokenize(lattice_view, tok_lattice);

  const char* prop_it = val_end;
  while (prop_it + prop_key_len <= end) {

    if (std::memcmp(prop_it, prop_key.data(), prop_key_len) == 0) {
      prop_it += prop_key_len;
      const char* prop_end = prop_it;
      while (prop_end < end && *prop_end != ' ')
        ++prop_end;
      std::string_view ppt = std::string_view(prop_it, prop_end);
      tokenizer_tmpl(ppt, tok_props, ColumnDelimiter{});
      break;
    }

    ++prop_it;
  }

  if (tok_props.size() == 0 || tok_props.size() % 3 != 0) return false;

  return true;
};

inline bool parse_exyz_lattice(const TokenSet& tokens, IOContext& ctx) {
  // read (ax ay az) (bx by bz) (cx cy cz)
  return tokens_to_num(
    tokens.at<0>(), ctx.H.m11,
    tokens.at<1>(), ctx.H.m21,
    tokens.at<2>(), ctx.H.m31,
    tokens.at<3>(), ctx.H.m12,
    tokens.at<4>(), ctx.H.m22,
    tokens.at<5>(), ctx.H.m32,
    tokens.at<6>(), ctx.H.m13,
    tokens.at<7>(), ctx.H.m23,
    tokens.at<8>(), ctx.H.m33
  );
}

template <size_t N> inline bool parse_exyz_properties(const TokenSetTmpl<N>& tokens, Properties& properties) {
  size_t begin = 0;
  size_t size = 0;
  for (size_t i = 0; i < tokens.size(); i += 3) {
    size_t j = i + 1;
    size_t k = i + 2;

    tokens_to_num(tokens[k], size);

    size_t index_field;
    size_t index_type;

    if (match_token_any(tokens[i], tok_props, index_field) && match_token_any(tokens[j], tok_field_type, index_type)) {
      properties.data[index_field] = {
          .field = props_list[index_field], .type = type_list[index_type], .size = size, .index = begin};
    }
    // TODO: validate type and size
    begin += size;
  }

  return true;
}

inline bool parse_line(IOContext& ctx, const Properties& properties, const TokenSet& tokens, size_t particle_count) {

  std::string_view type_as_string;
  double x;
  double y;
  double z;
  size_t type;
  size_t id;

  // parse type
  if (!ctx.species.get_type_from_sv(tokens[properties.species().index], type)) {
    lerr << "Undefined type: " << type_as_string << std::endl;
    return false;
  }

  // parse positions
  size_t index = properties.positions().index;
  if (!(tokens_to_num(tokens[index], x, tokens[index + 1], y, tokens[index + 2], z))) {
    lerr << "Unable to parse positions..." << std::endl;
    return false;
  }

  // optional atom id
  (properties.has_id()) ? tokens_to_num(tokens[properties.id().index], id) : id = particle_count;

  // TODO: optional velocity

  ParticleTupleIO& p = ctx.particle_data[particle_count];
  p[field::rx] = x;
  p[field::ry] = y;
  p[field::rz] = z;
  p[field::vx] = 0.0;
  p[field::vy] = 0.0;
  p[field::vz] = 0.0;
  p[field::id] = id;
  p[field::type] = type;

  return true;
}

class ExtendedXYZParser : public TextParser {
public:
  ExtendedXYZParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override;
  bool parse(IOContext& ctx) override;

private:
  TokenSet64 m_ptokens;
  Properties m_properties;
};

}

template<> inline const Metadata get_parser_metadata<xyz::ExtendedXYZParser>() {
  return {
    .name        = "Extended XYZ",
    .aliases     = { "xyz" },
    .description = "Extended XYZ"
  };
}

/* ------------------------------------------------------------------------- */

struct FileInfos {
  FileCompression compression = FileCompression::NONE;
  const RegisteredFormat format{};
};

std::string remove_leading_chars(const std::string& str, const char c);
const FileInfos guess_file_infos(const std::string& filepath);
const FileInfos get_file_infos(const std::string& filepath, std::string format, std::string compression);

} // namespace ioextra

} // namespace exaStamp
