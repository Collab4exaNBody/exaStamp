#pragma once

#include <array>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <utility>

#include <onika/file_utils.h>
#include <onika/log.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/string_utils.h>

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>

#include <exaStamp/io/tokenizer_utils.h>

#define READ_BUFFER_SIZE 1'048'576
#define FILE_BUFFER_SIZE READ_BUFFER_SIZE
#define LZMA_STREAM_BUF_SIZE READ_BUFFER_SIZE
#define LZMA_INTERNAL_BUF_SIZE READ_BUFFER_SIZE
#define BZ2_STREAM_BUF_SIZE READ_BUFFER_SIZE
#define BZ2_INTERNAL_BUF_SIZE READ_BUFFER_SIZE
#define NULL_CHAR '\0'

#define PANIC(...)                                                                                                     \
  do {                                                                                                                 \
    onika::lout << "ERR: [" << std::filesystem::path(__FILE__).filename().c_str() << ":" << __func__ << ":"            \
                << __LINE__ << "] " << onika::format_string(__VA_ARGS__) << std::endl;                                 \
    std::abort();                                                                                                      \
  } while (0);

#define linfo(...) onika::lout << onika::format_string(__VA_ARGS__);

namespace exaStamp {
using namespace exanb;

using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_vx, field::_vy, field::_vz,
                                                 field::_id, field::_type>;
using ParticleData = std::vector<ParticleTupleIO>;

namespace matrix_traits {

struct Mat3dSlice {
  double &x, &y, &z;
};

template <size_t I, size_t J>
constexpr inline double& at(Mat3d& m) {
  double* base = &m.m11;
  static_assert(I >= 0 && I < 3);
  static_assert(J >= 0 && J < 3);
  return base[I * 3 + J];
}

template <size_t N> inline constexpr Mat3dSlice column(Mat3d& m) {
  static_assert(N >= 0 && N < 2);
  double* base = &m.m11;
  return Mat3dSlice{base[N], base[N + 3], base[N + 6]};
}

template <size_t N> inline constexpr Mat3dSlice row(Mat3d& m) {
  static_assert(N >= 0 && N < 2);
  double* base = &m.m11;
  size_t index = 3 * N;
  return Mat3dSlice{base[index], base[index + 1], base[index + 2]};
}

inline Mat3dSlice column(size_t i, Mat3d& m) {
  double* base = &m.m11;
  return Mat3dSlice{base[i], base[i + 3], base[i + 6]};
}

inline Mat3dSlice row(size_t i, Mat3d& m) {
  double* base = &m.m11;
  size_t index = 3 * i;
  return Mat3dSlice{base[index], base[index + 1], base[index + 2]};
}

template <typename T>
concept VecLike = requires(T v) {
  { v.x } -> std::convertible_to<double>;
  { v.y } -> std::convertible_to<double>;
  { v.z } -> std::convertible_to<double>;
};
} // namespace matrix_traits

namespace ioextra {

using namespace tokens;
using namespace tokens::numeric;
using namespace enum_traits;

/* ------------------------------------------------------------------------- */

struct SpeciesMap {

  static inline uint32_t pack_species(StringView sv) {
    char a = sv.size() > 0 ? static_cast<uint8_t>(sv[0]) : 0;
    char b = sv.size() > 1 ? static_cast<uint8_t>(sv[1]) : 0;
    char c = sv.size() > 2 ? static_cast<uint8_t>(sv[2]) : 0;
    char d = sv.size() > 3 ? static_cast<uint8_t>(sv[3]) : 0;
    return (static_cast<uint32_t>(a) << 24 | static_cast<uint32_t>(b) << 16 | static_cast<uint32_t>(c) << 8 |
            static_cast<uint32_t>(d));
  }

  void add(const std::string& symbol, ssize_t index) {
    if (count >= 10)
      return;
    uint32_t key = pack_species(symbol);
    keys[count] = key;
    values[count] = index;
    count++;
  }

  ssize_t lookup(StringView symbol) {
    uint32_t key = pack_species(symbol);

    if (key == cache_key)
      return cache_index;

    ssize_t index = -1;
    for (size_t i = 0; i < count; ++i) {
      if (keys[i] == key) {
        index = values[i];
        break;
      }
    }

    if (cmp::lt(index, -1)) {
      cache_key = key;
      cache_index = index;
    }

    return index;
  }

  size_t count = 0;

  std::array<uint32_t, 10> keys{};
  std::array<ssize_t, 10> values{};
  uint32_t cache_key = 0xffffffff;
  int cache_index = -1;
};

struct IOContext {

  struct DataBuffer {
    size_t nat;
    size_t ntype;

    Mat3d cell{};
    Vec3d origin{};

    SpeciesMap species;
    ParticleData particles;
  };

  enum ContextFlags : uint16_t {
    REMAP_ATOM = 1 << 1,
    DOMAIN_ONLY = 1 << 2,
    TRICLINIC = 1 << 3,
  };

  uint16_t flags = 0;
  DataBuffer data{};
  Vec3d timer{};

  template <uint16_t Flag>
  inline bool has() {
    return (flags & (Flag));
  }

  template<int16_t Flag>
  void set(bool toggle) {
    (toggle) ? flags |= Flag : flags &= Flag;
  }

};

/* ------------------------------------------------------------------------- */

// File openning mode
enum FileMode : char {
  READ = 'r',
  WRITE = 'w',
  APPEND = 'a',
};

// File compression mode
enum FileCompression {
  NONE,
  GZIP,
  BZIP2,
  XZ,
};

/* ------------------------------------------------------------------------- */

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
    PANIC("Invalid convertion from Filemode to char *");
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

typedef struct gzFile_s* gzFile;

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
      : m_path(std::move(path)),
        m_mode(mode),
        m_compression(compression) {}

public:
  inline const std::string& path() const { return m_path; };
  inline FileMode mode() const { return m_mode; };
  inline FileCompression compression() const { return m_compression; };
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

  inline bool is_buffer_init() const { return m_buf[0] != NULL_CHAR; }
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

template <typename T>
const Metadata get_parser_metadata() {
  return {};
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

  virtual inline bool operator()(IOContext&) { return false; };
  virtual inline bool operator()(IOContext&, size_t) { return false; };
  virtual size_t size() = 0;
};

class TextParser : public Parser {

protected:
  bool m_eof = false;
  std::vector<size_t> m_loc{};
  TextFile m_file;

  std::array<std::string_view, 10> m_line_buffer{};
  TokenSet m_tokens{};

  inline std::string_view& current_line(void) { return m_line_buffer[0]; };
  inline const std::string_view& current_line(void) const { return m_line_buffer[0]; }

public:
  TextParser(std::string filepath, FileMode mode, FileCompression compression)
      : m_file(std::move(filepath), mode, compression) {}

  virtual ~TextParser() override = default;

  virtual int64_t next() = 0;
  virtual bool parse(IOContext& ctx) = 0;

  bool operator()(IOContext& ctx) override;
  bool operator()(IOContext& ctx, size_t index) override;

  void scan();

  inline bool eof() { return m_eof; }
  inline size_t size() override { return m_loc.size(); };
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

  template <typename T>
  void register_format() {

    size_t index = index_to_parser.size();
    Metadata metadata = get_parser_metadata<T>();

    index_to_parser.insert(
        {index,
         {.metadata = metadata,
          .creator = [](std::string filepath, FileMode mode, FileCompression compression) {
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
      return {.metadata = Metadata{}, .creator = nullptr, .is_valid = false};
    }
    return index_to_parser.find(it->second)->second;
  }

  inline bool has_alias(const std::string& alias) const {
    return alias_to_index.find(alias) != alias_to_index.end();
  }

  std::string infos() const;
};

const ParserFactory& parser_factory();

/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Data format

namespace lammps_data {

using AtomDataParserFn = bool (*)(const TokenSet&, IOContext&, size_t);
using KeywordParserFn = bool (*)(const TokenSet&, IOContext&);

enum class AtomStyle {
  Atomic,
  Charge,
  Molecular,
  Full,
  Count,
  None,
};

constexpr TokenNeedles<enum_size<AtomStyle>()> NEEDLES_ATOM_STYLE{
    "atomic",
    "charge",
    "molecule",
    "full",
};

enum class Keyword {
  Atoms,
  AtomTypes,
  XAxis,
  YAxis,
  ZAxis,
  Tilts,
  Section,
  Count,
};

// Needles for header keywords (HDRKEY)
constexpr TokenNeedles<1> NEEDLES_KEYWORD_ATOMS{"atoms"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_ATOM_TYPES{"atom", "types"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_XAXIS{"xlo", "xhi"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_YAXIS{"ylo", "yhi"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_ZAXIS{"zlo", "zhi"};
constexpr TokenNeedles<3> NEEDLES_KEYWORD_TILTS{"xy", "xz", "yz"};

enum class Section {
  Masses,
  Atoms,
  Velocities,
  Bonds,
  Angles,
  Dihedrals,
  Count,
  Header,
  Ignored,
};

constexpr TokenNeedles<enum_size<Section>()> NEEDLES_SECTIONS{
    "Masses",
    "Atoms",
    "Velocities",
    "Bonds",
    "Angles",
    "Dihedrals",
};

bool parse_keyword_atoms(const TokenSet&, IOContext&);
bool parse_keyword_atom_types(const TokenSet&, IOContext&);
bool parse_keyword_xaxis(const TokenSet&, IOContext&);
bool parse_keyword_yaxis(const TokenSet&, IOContext&);
bool parse_keyword_zaxis(const TokenSet&, IOContext&);
bool parse_keyword_tilts(const TokenSet&, IOContext&);
bool parse_keyword_section(const TokenSet&, IOContext&);

constexpr std::array<KeywordParserFn, enum_size<Keyword>()> keyword_parsers{
    &parse_keyword_atoms,
    &parse_keyword_atom_types,
    &parse_keyword_xaxis,
    &parse_keyword_yaxis,
    &parse_keyword_zaxis,
    &parse_keyword_tilts,
    &parse_keyword_section,
};

template <bool RemapAtomID>
inline bool parse_atom_data_atomic(const TokenSet& tokens, IOContext& ctx, size_t particle_count) {

  size_t index = 0;
  size_t id, type;
  double x, y, z;

  if (!tokens.atleast(5)) // atom-id atom-type x y z
    return false;
  if (!parse_single_token(tokens.at<0>(), id))
    return false;
  if (!parse_single_token(tokens.at<1>(), type))
    return false;
  if (!parse_single_token(tokens.at<2>(), x))
    return false;
  if (!parse_single_token(tokens.at<3>(), y))
    return false;
  if (!parse_single_token(tokens.at<4>(), z))
    return false;

  if constexpr (RemapAtomID) {
    index = particle_count - 1;
  } else {
    if (cmp::le(id, ctx.data.nat)) {
      index = id - 1;
    } else {
      PANIC("ID=%ld is out-of-range N=%ld", id, ctx.data.nat);
    }
  }

  ParticleTupleIO& p = ctx.data.particles[index];
  p[field::rx] = x;
  p[field::ry] = y;
  p[field::rz] = z;
  p[field::vx] = 0.0;
  p[field::vy] = 0.0;
  p[field::vz] = 0.0;
  p[field::id] = index;
  p[field::type] = type - 1;

  return true;
}

template <bool RemapAtomID>
inline bool parse_atom_data_charge(const TokenSet& tokens, IOContext& ctx, size_t particle_count) {

  size_t index = 0;
  size_t id, type;
  double x, y, z;

  if (!tokens.atleast(6)) // atom-id atom-type q x y z
    return false;
  if (!parse_single_token(tokens.at<0>(), id))
    return false;
  if (!parse_single_token(tokens.at<1>(), type))
    return false;
  if (!parse_single_token(tokens.at<3>(), x))
    return false;
  if (!parse_single_token(tokens.at<4>(), y))
    return false;
  if (!parse_single_token(tokens.at<5>(), z))
    return false;

  if constexpr (RemapAtomID) {
    index = particle_count - 1;
  } else {
    if (cmp::le(id, ctx.data.nat)) {
      index = id - 1;
    } else {
      PANIC("ID=%ld is out-of-range N=%ld", id, ctx.data.nat);
    }
  }

  ParticleTupleIO& p = ctx.data.particles[index];
  p[field::rx] = x;
  p[field::ry] = y;
  p[field::rz] = z;
  p[field::vx] = 0.0;
  p[field::vy] = 0.0;
  p[field::vz] = 0.0;
  p[field::id] = index;
  p[field::type] = type - 1;

  return true;
}

template <bool RemapAtomID>
inline bool parse_atom_data_molecular(const TokenSet& tokens, IOContext& ctx, size_t particle_count) {
  // [[maybe_unused]] constexpr size_t len = 6; // molecular: atom-id molecule-id atom-type x y z
  PANIC("Not Implemented style=molecular");
  return false;
}

template <bool RemapAtomID>
inline bool parse_atom_data_full(const TokenSet& tokens, IOContext& ctx, size_t particle_count) {
  // [[maybe_unused]] constexpr size_t len = 7; // full:      atom-id molecule-id atom-type q x y z
  PANIC("Not Implemented style=full");
  return false;
}

template <bool RemapAtomID>
inline const std::array<const AtomDataParserFn, enum_size<AtomStyle>()>& atom_data_parsers() {
  static constexpr std::array<const AtomDataParserFn, enum_size<AtomStyle>()> fptr_hdlr{
      &parse_atom_data_atomic<RemapAtomID>,
      &parse_atom_data_charge<RemapAtomID>,
      &parse_atom_data_molecular<RemapAtomID>,
      &parse_atom_data_full<RemapAtomID>,
  };
  return fptr_hdlr;
}

template <bool RemapAtomID>
inline bool parse_atom_velocities(const TokenSet& tokens, IOContext& ctx, size_t particle_count) {
  size_t id, index = 0;
  double vx, vy, vz;

  if (!tokens.atleast(4)) // atom-id vx vy vz
    return false;
  if (!parse_single_token(tokens.at<0>(), id))
    return false;
  if (!parse_single_token(tokens.at<3>(), vx))
    return false;
  if (!parse_single_token(tokens.at<4>(), vy))
    return false;
  if (!parse_single_token(tokens.at<5>(), vz))
    return false;

  if constexpr (RemapAtomID) {
    index = particle_count - 1;
  } else {
    if (cmp::le(id, ctx.data.nat)) {
      index = id - 1;
    } else {
      PANIC("ID=%ld is out-of-range N=%ld", id, ctx.data.nat);
    }
  }

  ParticleTupleIO& p = ctx.data.particles[index];
  p[field::vx] = vx;
  p[field::vy] = vy;
  p[field::vz] = vz;
  return true;
}

class LAMMPSDataParser : public TextParser {
public:
  LAMMPSDataParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression), m_atom_style(AtomStyle::None), m_current_section(Section::Header) {}

  int64_t next() override;
  bool parse(IOContext&) override;

private:
  AtomStyle m_atom_style;
  Section m_current_section;

  // keep track of the current section line
  inline std::string_view& section_line(void) { return m_line_buffer[1]; }
  inline const std::string_view& section_line(void) const { return m_line_buffer[1]; }

  bool parse_section_header(IOContext& ctx);
  bool parse_section_atoms(IOContext& ctx);
  bool parse_section_masses(void);
  bool parse_section_velocities(IOContext& ctx);
  void jumpto_next_section();

  inline void set_current_section(std::string_view str) {
    size_t index;
    if (token_match_any(str, NEEDLES_SECTIONS, index)) {
      m_current_section = enum_from_index<Section>(index);
    } else {
      m_current_section = Section::Ignored;
    }
  }
};

} // namespace lammps_data

template <>
inline const Metadata get_parser_metadata<lammps_data::LAMMPSDataParser>() {
  return {
      .name = "LAMMPS Data",
      .aliases = {"lmp", "lmp-data", "data"},
      .description = "LAMMPS Data format",
  };
}

/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Data format

namespace lammps_dump {

enum class Field {
  ID,      // atom id
  MolID,   // molecule id
  Type,    // atom type
  Element, // element
  X,       // x wrapped
  Y,       // y wrapped
  Z,       // z wrapped
  XU,      // x unwrapped
  YU,      // y unwrapped
  ZU,      // z unwrapped
  XS,      // x scaled position
  YS,      // y scaled position
  ZS,      // z scaled position
  XSU,     // x scaled unwrapped position
  YSU,     // y scaled unwrapped position
  ZSU,     // z scaled unwrapped position
  VX,      // x velocity
  VY,      // y velocity
  VZ,      // z velocity
  None,    // ignored field
  Count,   //
};

constexpr TokenNeedles<enum_size<Field>()> NEEDLES_FIELDS = {
    "id",
    "mol",
    "type",
    "element",
    "x",
    "y",
    "z",
    "xu",
    "yu",
    "zu",
    "xs",
    "ys",
    "zs",
    "xsu",
    "ysu",
    "zsu",
    "vx",
    "vy",
    "vz",
};

struct Coordinates {
  bool scaled;
  bool unwrapped;
  Field x, y, z;
};

constexpr std::array<Coordinates, 4> COORDINATES{{
    {false, false, Field::X, Field::Y, Field::Z},
    {false, true, Field::XU, Field::YU, Field::ZU},
    {true, false, Field::XS, Field::YS, Field::ZS},
    {true, true, Field::XSU, Field::YS, Field::ZS},
}};

struct Property {
  Field field = Field::None;
  size_t index = 0;
};

struct Properties : properties_traits::GenericContainer<Property, Field, &Property::field> {
  bool triclinic = false;
  bool scaled = false;
  bool unwrapped = false;

  IJK loc_positions{};
  IJK loc_velocities{};

  bool (*parse_id)(size_t&, const TokenSet&, size_t, size_t) = nullptr;
  bool (*parse_molid)(size_t&, const TokenSet&, size_t) = nullptr;
  bool (*parse_type)(size_t&, const TokenSet&, size_t) = nullptr;
  bool (*parse_positions)(Vec3d&, const TokenSet&, const IJK&, const IOContext&) = nullptr;
  bool (*parse_velocities)(Vec3d&, const TokenSet&, const IJK&) = nullptr;

  template <Field i, Field j, Field k> inline constexpr IJK ijk() const {
    return IJK{
        static_cast<ssize_t>(get<i>().index),
        static_cast<ssize_t>(get<j>().index),
        static_cast<ssize_t>(get<k>().index),
    };
  }

  inline IJK ijk(Field i, Field j, Field k) {
    return {
        static_cast<ssize_t>(this->get(i).index),
        static_cast<ssize_t>(this->get(j).index),
        static_cast<ssize_t>(this->get(k).index),
    };
  }

  inline bool has_positions() {
    return (has<Field::X, Field::Y, Field::Z>() || has<Field::XU, Field::YU, Field::ZU>() ||
            has<Field::XS, Field::YS, Field::ZS>() || has<Field::XSU, Field::YSU, Field::ZSU>());
  }
};

template <bool IsTriclinic, matrix_traits::VecLike T>
inline bool parse_box_bounds_axis(const TokenSet& tokens, const T& v) {
  if constexpr (IsTriclinic) {
    if (!tokens.atleast(3))
      return false;
    if (!parse_single_token(tokens.at<0>(), v.x))
      return false;
    if (!parse_single_token(tokens.at<1>(), v.y))
      return false;
    if (!parse_single_token(tokens.at<2>(), v.z))
      return false;
    return true;
  } else {
    if (!tokens.atleast(2))
      return false;
    if (!parse_single_token(tokens.at<0>(), v.x))
      return false;
    if (!parse_single_token(tokens.at<1>(), v.y))
      return false;
    return true;
  }
}

inline void convert_bbox_to_lattice_vectors(IOContext& ctx, const Mat3d& bbox) {
  double xlo = bbox.m11 - std::min(bbox.m13 + bbox.m23, std::min(0., bbox.m23));
  double xhi = bbox.m12 - std::max(bbox.m13 + bbox.m23, std::max(0., bbox.m23));
  double ylo = bbox.m21 - std::min(0., bbox.m33);
  double yhi = bbox.m22 - std::max(0., bbox.m33);
  double zlo = bbox.m31;
  double zhi = bbox.m32;

  ctx.data.cell.m11 = xhi - xlo;
  ctx.data.cell.m22 = yhi - ylo;
  ctx.data.cell.m33 = zhi - zlo;
  ctx.data.cell.m12 = bbox.m13;
  ctx.data.cell.m13 = bbox.m23;
  ctx.data.cell.m23 = bbox.m33;

  ctx.data.origin.x = xlo;
  ctx.data.origin.y = ylo;
  ctx.data.origin.z = zlo;
};

template <bool Scaled, bool Unwrapped, bool Triclinic>
inline bool parse_field_positions(Vec3d& r, const TokenSet& tokens, const IJK& inds, const IOContext& ctx) {
  if (!parse_single_token(tokens[inds.i], r.x))
    return false;
  if (!parse_single_token(tokens[inds.j], r.y))
    return false;
  if (!parse_single_token(tokens[inds.k], r.z))
    return false;

  if constexpr (Scaled) {
    if constexpr (Triclinic) {
      // x = xlo + x * lx + y * xy + z * xz
      // y = zlo + y * ly + z * xz
      // z = ylo + z * lz;
      r.x = ctx.data.origin.x + r.x * ctx.data.cell.m11 + r.y + ctx.data.cell.m12 + r.z * ctx.data.cell.m13;
      r.y = ctx.data.origin.y + r.y * ctx.data.cell.m22 + r.z * ctx.data.cell.m13;
      r.z = ctx.data.origin.z + r.z * ctx.data.cell.m33;
    } else {
      r.x = ctx.data.origin.x + r.x * ctx.data.cell.m11;
      r.y = ctx.data.origin.y + r.y * ctx.data.cell.m22;
      r.z = ctx.data.origin.z + r.z * ctx.data.cell.m33;
    }
  }
  return true;
}

template <bool HasID>
inline bool parse_field_id(size_t& id, const TokenSet& tokens, size_t index, size_t particle_count) {
  if constexpr (HasID) {
    return parse_single_token(tokens[index], id);
  } else {
    id = particle_count;
    return true;
  }
}

template <bool HasMolID>
inline bool parse_field_molid(size_t& molid, const TokenSet& tokens, size_t index) {
    if constexpr (HasMolID) {
        return parse_single_token(tokens[index], molid);
    } else {
        molid = 0;
        return true;
    }
}

template <bool HasParticleType>
inline bool parse_field_type(size_t& type, const TokenSet& tokens, size_t index) {
  if constexpr (HasParticleType) {
    return parse_single_token(tokens[index], type);
  } else {
    type = 0;
    return true;
  }
}

template <bool HasVelocities>
inline bool parse_field_velocities(Vec3d& v, const TokenSet& tokens, const IJK& inds) {
  if constexpr (HasVelocities) {
    if (!parse_single_token(tokens[inds.i], v.x))
      return false;
    if (!parse_single_token(tokens[inds.j], v.y))
      return false;
    if (!parse_single_token(tokens[inds.k], v.z))
      return false;
    return true;
  } else {
    v = {0., 0., 0};
    return true;
  }
}

inline bool parse_atom_properties(const TokenSet& tokens, Properties& properties) {

  constexpr size_t headlen = 2; // ITEM: ATOMS
  size_t ifield;

  for (size_t i = headlen; i < tokens.size(); ++i) {
    if (token_match_any(tokens[i], NEEDLES_FIELDS, ifield)) {
      properties.set({.field = enum_from_index<Field>(ifield), .index = i - headlen});
    }
  }

  properties.parse_id = (properties.has<Field::ID>()) ? parse_field_id<true> : parse_field_id<false>;
  properties.parse_molid = (properties.has<Field::MolID>()) ? parse_field_molid<true> : parse_field_molid<false>;
  properties.parse_type = (properties.has<Field::Type>()) ? parse_field_type<true> : parse_field_type<false>;

  // deduced coodinate system base on present field
  for (const Coordinates& c : COORDINATES) {
    if (properties.has(c.x, c.y, c.z)) {
      properties.scaled = c.scaled;
      properties.unwrapped = c.unwrapped;
      properties.loc_positions = properties.ijk(c.x, c.y, c.z);
      properties.expect(c.x, c.y, c.z);
    }
  }

  if (!properties.has_positions()) {
    PANIC("No positions fields are present...");
  }

  if (properties.scaled && properties.triclinic) {
    properties.parse_positions = parse_field_positions<true, false, true>;
  } else if (properties.scaled && !properties.triclinic) {
    properties.parse_positions = parse_field_positions<true, false, false>;
  } else if (!properties.scaled && properties.triclinic) {
    properties.parse_positions = parse_field_positions<false, false, true>;
  } else {
    properties.parse_positions = parse_field_positions<false, false, false>;
  }

  // check if there are velocity fields
  if (properties.has<Field::VX, Field::VY, Field::VZ>()) {
    properties.parse_velocities = parse_field_velocities<true>;
    properties.loc_velocities = properties.ijk<Field::VX, Field::VY, Field::VZ>();
    properties.expect(Field::VX, Field::VY, Field::VZ);
  } else {
    properties.parse_velocities = parse_field_velocities<false>;
  }

  return true;
}

template <bool RemapAtomID>
inline bool parse_atom_data(const TokenSet& tokens, IOContext& ctx, const Properties& ppt, size_t particle_count) {

  size_t index = 0;
  size_t id, molid, type;
  Vec3d pos, vel;

  if (!tokens.atleast(ppt.len()))
    return false;
  if (!ppt.parse_id(id, tokens, ppt.get<Field::ID>().index, particle_count))
    return false;
  if (!ppt.parse_molid(molid, tokens, ppt.get<Field::MolID>().index))
    return false;
  if (!ppt.parse_type(type, tokens, ppt.get<Field::Type>().index))
    return false;
  if (!ppt.parse_positions(pos, tokens, ppt.loc_positions, ctx))
    return false;
  if (!ppt.parse_velocities(vel, tokens, ppt.loc_velocities))
    return false;

  if constexpr (RemapAtomID) {
    index = particle_count - 1;
  } else {
    if (cmp::le(id, ctx.data.nat)) {
      index = id - 1;
    } else {
      PANIC("ID=%ld is out-of-range N=%ld", id, ctx.data.nat);
    }
  }

  ParticleTupleIO& p = ctx.data.particles[index];
  p[field::rx] = pos.x;
  p[field::ry] = pos.y;
  p[field::rz] = pos.z;
  p[field::vx] = vel.x;
  p[field::vy] = vel.y;
  p[field::vz] = vel.z;
  p[field::id] = index;
  p[field::type] = type - 1;

  return true;
}

class LAMMPSDumpParser : public TextParser {
public:
  LAMMPSDumpParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override;
  bool parse(IOContext& ctx) override;
};

} // namespace lammps_dump

template <>
inline const Metadata get_parser_metadata<lammps_dump::LAMMPSDumpParser>() {
  return {
      .name = "LAMMPS Dump",
      .aliases = {"dump", "lmp-dump", "lammps-dump"},
      .description = "LAMMPS Dump format",
  };
}

namespace xyz {

enum class FieldType : uint8_t {
  R = 1 << 0,
  S = 1 << 1,
  I = 1 << 2,
  Count = 3,
  None = 0,
};

constexpr TokenNeedles<enum_size<FieldType>()> NEEDLES_FIELD_TYPE{
    "R",
    "S",
    "I",
};

enum class Field {
  Species,
  Positions,
  Velocities,
  ID,
  Type,
  None,
  Count,
};

constexpr TokenNeedles<enum_size<Field>()> NEEDLES_FIELDS{
    "species",
    "pos",
    "velo",
    "id",
    "type",
};

struct Property {
  Field field = Field::None;
  FieldType type = FieldType::None;
  size_t icol = 0;
  size_t ncol = 0;
};

struct Properties : public properties_traits::GenericContainer<Property, Field, &Property::field> {

  using Enum = typename properties_traits::GenericContainer<Property, Field, &Property::field>::enum_t;

  template <Enum... args> constexpr void expect() {
    // static_assert((std::is_same_v<Enum, Enum> && args...), "Only E enum type allowed");
    ((length += (has<args>()) ? get<args>().ncol : static_cast<size_t>(0)), ...);
  }

  IJK loc_positions{};
  IJK loc_velocities{};

  bool (*parse_id)(size_t&, const TokenSet&, size_t, size_t) = nullptr;
  bool (*parse_type)(ssize_t&, const TokenSet&, size_t, size_t, IOContext&) = nullptr;
  bool (*parse_velocities)(Vec3d&, const TokenSet&, const IJK&) = nullptr;
};

template <size_t N1, size_t N2>
bool parse_extxyz_comment_line(StringView line, TokenSetTmpl<N1>& lattice_tokens, TokenSetTmpl<N2>& properties_tokens) {
  constexpr StringView key_lattice = "Lattice=\"";
  constexpr size_t key_lattice_len = key_lattice.size(); // 9
  constexpr StringView key_properties = "Properties=";
  constexpr size_t key_properties_len = key_properties.size(); // 11

  const char* begin = line.data();
  const char* end = line.data() + line.size();

  if (!startswith(line, key_lattice)) {
    return false;
  }

  // 1. Parse Lattice=...
  const char* first = begin + key_lattice_len;
  const char* last = first;

  while (last < end && !DoubleQuoteDelimiter{}(last)) {
    ++last;
  }

  if (last >= end) {
    return false;
  }

  tokenize(StringView(first, last), lattice_tokens);

  // 2. Parse Properties=....
  first = last;
  while (first + key_properties_len <= end) {

    if (std::memcmp(first, key_properties.data(), key_properties_len) == 0) {
      first += key_properties_len;
      last = first;

      while (last < end && !SpaceDelimiter{}(last)) {
        ++last;
      }

      tokenize_impl(StringView(first, last), properties_tokens, ColumnDelimiter{});
      break;
    }
    ++first;
  }

  return true;
}

inline bool parse_extxyz_lattice(const TokenSet& tokens, IOContext& ctx) {
  if (!tokens.atleast(9))
    return false;
  if (!parse_single_token(tokens.at<0>(), ctx.data.cell.m11)) // ax
    return false;
  if (!parse_single_token(tokens.at<1>(), ctx.data.cell.m21)) // ay
    return false;
  if (!parse_single_token(tokens.at<2>(), ctx.data.cell.m31)) // az
    return false;
  if (!parse_single_token(tokens.at<3>(), ctx.data.cell.m12)) // bx
    return false;
  if (!parse_single_token(tokens.at<4>(), ctx.data.cell.m22)) // by
    return false;
  if (!parse_single_token(tokens.at<5>(), ctx.data.cell.m32)) // bz
    return false;
  if (!parse_single_token(tokens.at<6>(), ctx.data.cell.m13)) // cx
    return false;
  if (!parse_single_token(tokens.at<7>(), ctx.data.cell.m23)) // cy
    return false;
  if (!parse_single_token(tokens.at<8>(), ctx.data.cell.m33)) // cz
    return false;
  return true;
}

template <bool HasID>
inline bool parse_field_id(size_t& id, const TokenSet& tokens, size_t index, size_t particle_count) {
  if constexpr (HasID) {
    return parse_single_token(tokens[index], id);
  } else {
    id = particle_count;
    return true;
  }
}

template <bool HasMolID>
inline bool parse_field_molid(size_t& molid, const TokenSet& tokens, size_t index) {
    if constexpr (HasMolID) {
        return parse_single_token(tokens[index], molid);
    } else {
        molid = 0;
        return true;
    }
}


template <bool HasParticleSpecies, bool HasParticleType>
inline bool parse_field_type(ssize_t& type, const TokenSet& tokens, size_t i, size_t j, IOContext& ctx) {
  if constexpr (HasParticleSpecies) {
    type = ctx.data.species.lookup(tokens[i]);
    return cmp::ne(type, -1);
  } else if constexpr (HasParticleType) {
    return parse_single_token(tokens[j], type);
  } else {
    type = 0;
    return true;
  }
}

template <bool HasVelocities>
inline bool parse_field_velocities(Vec3d& v, const TokenSet& tokens, const IJK& inds) {
  if constexpr (HasVelocities) {
    if (!parse_single_token(tokens[inds.i], v.x))
      return false;
    if (!parse_single_token(tokens[inds.j], v.y))
      return false;
    if (!parse_single_token(tokens[inds.k], v.z))
      return false;
    return true;
  } else {
    v = {0., 0., 0};
    return true;
  }
}

inline bool parse_field_positions(Vec3d& r, const TokenSet& tokens, const IJK& inds) {
  if (!parse_single_token(tokens[inds.i], r.x))
    return false;
  if (!parse_single_token(tokens[inds.j], r.y))
    return false;
  if (!parse_single_token(tokens[inds.k], r.z))
    return false;
  return true;
}

inline bool parse_xyz_properties(const TokenSet& tokens, Properties& properties) {

  if (tokens.size() % 3 != 0) {
    return false;
  }

  size_t icol = 0;
  size_t ncol, ifield, itype;

  properties.clear();

  for (size_t i = 0; i < tokens.size(); i += 3) {
    size_t j = i + 1;
    size_t k = i + 2;

    parse_single_token(tokens[k], ncol);

    if ((token_match_any(tokens[i], NEEDLES_FIELDS, ifield) && token_match_any(tokens[j], NEEDLES_FIELD_TYPE, itype))) {
      properties.set({
          .field = enum_from_index<Field>(ifield),
          .type = enum_from_index<FieldType>(itype),
          .icol = icol,
          .ncol = ncol,
      });
    }

    icol += ncol;
  }

  // Field validation
  if (properties.has<Field::ID>()) {
    properties.parse_id = parse_field_id<true>;
    properties.expect<Field::ID>();
  } else {
    properties.parse_id = parse_field_id<false>;
  }

  if (!(properties.has<Field::Type>() || properties.has<Field::Species>())) {
    PANIC("Either Species or TYPE should be define in xyz properties")
  } else if (properties.has<Field::Species>()) {
    properties.parse_type = parse_field_type<true, false>;
    properties.expect<Field::Species>();
  } else if (properties.has<Field::Type>()) {
    properties.parse_type = parse_field_type<false, true>;
    properties.expect<Field::Type>();
  }

  if (!properties.has<Field::Positions>()) {
    PANIC("Positions are not defined in atom properties");
  }

  ssize_t ipos = static_cast<ssize_t>(properties.get<Field::Positions>().icol);
  properties.loc_positions = {ipos, ipos + 1, ipos + 2};

  if (properties.has<Field::Velocities>()) {
    properties.parse_velocities = parse_field_velocities<true>;
    ssize_t ivel = static_cast<ssize_t>(properties.get<Field::Velocities>().icol);
    properties.loc_velocities = {ivel, ivel + 1, ivel + 2};
  } else {
    properties.parse_velocities = parse_field_velocities<false>;
  }
  return true;
}

template <bool RemapAtomID>
inline bool parse_extxyz_atom_line(const TokenSet& tokens, IOContext& ctx, const Properties& ppt,
                                   size_t particle_count) {

  size_t index = 0;
  size_t id;
  ssize_t type;
  Vec3d pos, vel;

  if (!tokens.atleast(ppt.len())) {
    linfo("ERR: line too shorts %ld %ld\n", tokens.size(), ppt.len());
    return false;
  }
  if (!parse_field_positions(pos, tokens, ppt.loc_positions)) {
    linfo("ERR: failed to parse positions\n");
    return false;
  }
  if (!ppt.parse_id(id, tokens, ppt.get<Field::ID>().icol, particle_count)) {
    linfo("ERR: failed to parse atom id\n");
    return false;
  }
  if (!ppt.parse_type(type, tokens, ppt.get<Field::Species>().icol, ppt.get<Field::Type>().icol, ctx)) {
    linfo("ERR: failed to parse atom type\n");
    return false;
  }
  if (!ppt.parse_velocities(vel, tokens, ppt.loc_velocities)) {
    linfo("ERR: failed to parse velocities\n");
    return false;
  }

  if constexpr (RemapAtomID) {
    index = particle_count - 1;
  } else {
    if (cmp::le(id, ctx.data.nat)) {
      index = id - 1;
    } else {
      PANIC("ID=%ld is out-of-range N=%ld", id, ctx.data.nat);
    }
  }

  ParticleTupleIO& p = ctx.data.particles[index];
  p[field::rx] = pos.x;
  p[field::ry] = pos.y;
  p[field::rz] = pos.z;
  p[field::vx] = vel.x;
  p[field::vy] = vel.y;
  p[field::vz] = vel.z;
  p[field::id] = index;
  p[field::type] = static_cast<size_t>(type - 1);

  return true;
}

class ExtendedXYZParser : public TextParser {
public:
  ExtendedXYZParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override;
  bool parse(IOContext& ctx) override;

private:
  TokenSet m_ptokens;
  Properties m_properties;
};

} // namespace xyz

template <>
inline const Metadata get_parser_metadata<xyz::ExtendedXYZParser>() {
  return {
    .name        = "Extended XYZ",
    .aliases     = { "xyz", "exyz", "ext-xyz" },
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
bool read_external_atom_file(IOContext& ctx, std::string filepath, std::string format, std::string compression);

} // namespace ioextra
} // namespace exaStamp

/* ------------------------------------------------------------------------- */

// // template <size_t N> struct SpeciesTmpl : TokenSetTmpl<N> {
// //   inline bool get_type_from_sv(std::string_view type, size_t& index) {
// //     return match_token_any(type, this->tokens, index);
// //   }
// // };

// // using Species = SpeciesTmpl<16>;

// // // ExaStamp context when parsing a file.
// // struct Context {

// //   ParticleData& particle_data;

// //   Species species;

// //   size_t n_particles = 0;
// //   size_t n_species = 0; // number of species detected in the file

// //   Mat3d H = make_identity_matrix(); // lattice vector
// //   Vec3d origin{0., 0., 0.};         // box origin

// //   Mat3d Hinv{};
// //   Mat3d D{};
// //   Mat3d D2{};
// //   Mat3d Dinv{};
// //   
// //   Mat3d Xform = Mat3d{};            // initial xform (likely == I)
// //   Mat3d invXform = Mat3d{};         // intial xform^-1
// //   Vec3d boxlen{};                   // box len in x, y, z

// //   bool uniform_scale;               //

// //   // for timing purpose
// //   size_t elapsed_time_ms = 0;
// //   size_t elapsed_time_ys = 0;
// //   size_t elapsed_time_ns = 0;

// //   inline void init_domain(Domain& domain, ReadBoundsSelectionMode bounds_mode) {

// //     // assume H follow column vector convention
// //     Vec3d a{H.m11, H.m21, H.m31};
// //     Vec3d b{H.m12, H.m22, H.m32};
// //     Vec3d c{H.m13, H.m23, H.m33};

// //     boxlen.x = norm(a);
// //     boxlen.y = norm(b);
// //     boxlen.z = norm(c);

// //     lout << "a = " << a << onika::format_string(" , norm = %10.5f ", boxlen.x) << std::endl;
// //     lout << "b = " << b << onika::format_string(" , norm = %10.5f ", boxlen.y) << std::endl;
// //     lout << "c = " << c << onika::format_string(" , norm = %10.5f ", boxlen.z) << std::endl;

// //     D = diag_matrix(boxlen);
// //     Dinv = inverse(D);
// //     Mat3d F1 = H * Dinv;

// //     double cellsize = domain.cell_size();

// //     // check for grid dims. If grid dims is not provided by the user, 
// //     // it will be deducde from box lengths and cellsize.
// //     IJK grid_dim = domain.grid_dimension();
// //     size_t nx = static_cast<ssize_t>(boxlen.x / cellsize);
// //     size_t ny = static_cast<ssize_t>(boxlen.y / cellsize);
// //     size_t nz = static_cast<ssize_t>(boxlen.z / cellsize);
// //     nx = (grid_dim.i > 0 && cmp::ne(grid_dim.i, nx)) ? grid_dim.i : nx;
// //     ny = (grid_dim.j > 0 && cmp::ne(grid_dim.j, ny)) ? grid_dim.j : ny;
// //     nz = (grid_dim.k > 0 && cmp::ne(grid_dim.k, nz)) ? grid_dim.k : nz;

// //     D2 = diag_matrix(Vec3d{nx * cellsize, ny * cellsize, nz * cellsize});
// //     Mat3d D2inv = inverse(D2);
// //     Mat3d F2 = D * D2inv;

// //     Xform = F1 * F2;
// //     invXform = inverse(Xform);

// //     AABB domain_bounds = {{0., 0., 0.}, {D2.m11, D2.m22, D2.m33}};
// //     compute_domain_bounds(domain, bounds_mode, 0.0, domain_bounds, domain_bounds, true);
// //     uniform_scale = is_uniform_scale(Xform);

// //     if (uniform_scale) {
// //       domain.set_xform(make_identity_matrix());
// //       domain.set_cell_size(domain.cell_size() * Xform.m11);
// //       domain.set_bounds({domain.origin() * Xform.m11, domain.extent() * Xform.m11});
// //     } else {
// //       domain.set_xform(Xform);
// //     }
// //   }
// // };

// // using IOContext = Context;

// // inline void wrap_to_domain(Vec3d& r, const IOContext& ctx) {
// //     r.x = wrap(r.x, ctx.D2.m11);
// //     r.y = wrap(r.y, ctx.D2.m22);
// //     r.z = wrap(r.z, ctx.D2.m33);
// // }

// struct IOContext {};

// /* ------------------------------------------------------------------------- */

// // File openning mode
// enum FileMode : char {
//   READ = 'r',
//   WRITE = 'w',
//   APPEND = 'a',
// };

// // File compression mode
// enum FileCompression {
//   NONE,
//   GZIP,
//   BZIP2,
//   XZ,
// };

// /* ------------------------------------------------------------------------- */

// inline const char* convert_mode_to_char(FileMode mode) {
//   switch (mode) {
//   case FileMode::READ:
//     return "rb";
//     break;
//   case FileMode::APPEND:
//     return "a+b";
//     break;
//   case FileMode::WRITE:
//     return "wb";
//     break;
//   default:
//     PANIC("Invalid convertion from Filemode to char *");
//     break;
//   }
// }

// inline FileCompression convert_char_to_compression(const std::string& str) {
//   if (str == "gz")
//     return FileCompression::GZIP;
//   if (str == "bz2")
//     return FileCompression::BZIP2;
//   if (str == "xz")
//     return FileCompression::XZ;
//   return FileCompression::NONE;
// }

// inline std::string convert_compression_to_char(FileCompression& compression) {
//   if (compression == FileCompression::GZIP)
//     return "gz";
//   if (compression == FileCompression::BZIP2)
//     return "bz2";
//   if (compression == FileCompression::XZ)
//     return "xz";
//   return "";
// }

// /* ------------------------------------------------------------------------- */

// // Base class that handler basic file operation
// // Specialization allow abstraction on file types.
// class TextFileHandler {
// public:
//   TextFileHandler(const std::string& filepath);

//   virtual ~TextFileHandler() = default;

//   virtual void clear() noexcept = 0;
//   virtual void seek(uint64_t position) = 0;
//   virtual size_t read(char* buffer, size_t size) = 0;

// protected:
//   const std::string& path() { return m_path; };
//   std::string m_path;
// };

// class ASCIIFileHandler final : public TextFileHandler {
// public:
//   ASCIIFileHandler(const std::string& path, FileMode mode);
//   ~ASCIIFileHandler() override;

//   void clear() noexcept override;
//   void seek(uint64_t cursor) override;

//   size_t read(char* data, size_t size) override;

// private:
//   std::FILE* m_fptr;
// };

// #ifdef USE_ZLIB
// #define ZLIB_CONST
// #include <zconf.h>
// #include <zlib.h>

// typedef struct gzFile_s* gzFile;

// class GzipFileHandler final : public TextFileHandler {
// public:
//   GzipFileHandler(const std::string& path, FileMode mode);
//   ~GzipFileHandler() override;

//   void clear() noexcept override;
//   void seek(uint64_t cursor) override;

//   size_t read(char* buffer, size_t size) override;

//   const char* gz_error() const;
//   static unsigned safe_cast(size_t);

// private:
//   gzFile m_fptr = nullptr;
// };

// #endif

// #ifdef USE_LZMA
// #include <lzma.h>

// class XzFileHandler final : public TextFileHandler {
// public:
//   XzFileHandler(const std::string& path, FileMode mode);
//   ~XzFileHandler() override;

//   void clear() noexcept override;
//   void seek(uint64_t cursor) override;

//   size_t read(char* buffer, size_t size) override;

//   size_t safe_cast(uint64_t value);
//   void check_lzma_ret(lzma_ret ret);
//   void start_lzma_decoder_stream(lzma_stream* stream);

// private:
//   std::FILE* m_fptr = nullptr;
//   lzma_stream m_xz_stream = LZMA_STREAM_INIT;

//   std::vector<uint8_t> m_xz_buffer = {0};
//   // static constexpr size_t m_xz_buffer_end = IOEXTRA_LZMA_INTERNAL_BUF_SIZE;
//   // uint8_t* m_xz_buffer[m_xz_buffer_end];
// };

// #endif

// #ifdef USE_BZIP2
// #include <bzlib.h>

// class Bzip2FileHandler final : public TextFileHandler {
// public:
//   Bzip2FileHandler(const std::string& path, FileMode mode);
//   ~Bzip2FileHandler() override;

//   void clear() noexcept override;
//   void seek(uint64_t cursor) override;

//   size_t read(char* buffer, size_t size) override;

//   unsigned safe_cast(uint64_t size);
//   void check_bz2_retcode(int code);

// private:
//   std::FILE* m_fptr;
//   std::function<int(bz_stream*)> m_end_bz2_stream;
//   bz_stream m_bz2_stream;
//   std::vector<char> m_bz2_buffer;
// };

// #endif

// /* ------------------------------------------------------------------------- */

// class File {
// protected:
//   std::string m_path;
//   FileMode m_mode;
//   FileCompression m_compression;

//   File(std::string path, FileMode mode, FileCompression compression)
//       : m_path(std::move(path)), m_mode(mode), m_compression(compression) {}

// public:
//   inline const std::string& path() const { return m_path; };
//   inline const FileMode mode() const { return m_mode; };
//   inline const FileCompression compression() const { return m_compression; };
// };

// class TextFile : public File {
// private:
//   std::vector<char> m_buf;
//   const char* m_lstart; // line start
//   const char* m_bend;   // buf end

//   uint64_t m_cursor = 0;
//   bool m_eof = false;
//   bool m_hdlr_eof = false;

//   std::unique_ptr<TextFileHandler> m_handler;

//   inline bool is_buffer_init() const { return m_buf[0] != NULL_CHAR; }
//   void fill_buffer(size_t pos);

// public:
//   TextFile(std::string filepath, FileMode mode, FileCompression compression);

//   inline bool eof() const { return m_eof; };

//   uint64_t tell() const;
//   void seek(uint64_t pos);
//   void clear();
//   void reset();

//   // get the next line in the buffer.
//   // If the line is not complete, load the next section of the buffer
//   std::string_view get_line();
// };

// /* ------------------------------------------------------------------------- */

// struct Metadata {
//   const std::string name = "";
//   const std::vector<std::string> aliases = {};
//   const std::string description = "";
// };

// template <typename T> const Metadata get_parser_metadata() { return {}; }

// /* ------------------------------------------------------------------------- */

// // Base class to all format parsers
// class Parser {
// public:
//   Parser() = default;
//   virtual ~Parser() = default;

//   Parser(const Parser&) = delete;
//   Parser& operator=(const Parser&) = delete;
//   Parser(Parser&&) = delete;
//   Parser& operator=(Parser&&) = delete;

//   // read the first step in the file
//   virtual inline bool operator()(IOContext& ctx) { return false; };
//   // read the step at the given index
//   virtual inline bool operator()(IOContext& ctx, size_t index) { return false; };
//   // returns the number of known step in the file
//   virtual size_t size() = 0;
// };

// // Specialization for format that are just text files
// class TextParser : public Parser {

// protected:
//   bool m_eof = false;
//   std::vector<size_t> m_loc{};
//   TextFile m_file;

//   std::array<std::string_view, 10> m_line_buffer{};
//   TokenSet m_tokens{};

//   inline std::string_view& current_line(void) { return m_line_buffer[0]; };
//   inline const std::string_view& current_line(void) const { return m_line_buffer[0]; }

// public:
//   TextParser(std::string filepath, FileMode mode, FileCompression compression)
//       : m_file(std::move(filepath), mode, compression) {}

//   virtual ~TextParser() override = default;

//   virtual int64_t next() = 0;
//   virtual bool parse(IOContext& ctx) = 0;

//   bool operator()(IOContext& ctx) override;
//   bool operator()(IOContext& ctx, size_t index) override;

//   void scan();

//   inline bool eof() { return m_eof; }
//   inline size_t size() override { return m_loc.size(); };
// };

// /* ------------------------------------------------------------------------- */

// using ParserCreator = std::function<std::unique_ptr<Parser>(std::string, FileMode, FileCompression)>;

// struct RegisteredFormat {
//   const Metadata metadata;
//   const ParserCreator creator;
//   bool is_valid = true;
// };

// class ParserFactory {
// public:
//   ParserFactory();

//   std::unordered_map<std::string, size_t> alias_to_index;
//   std::unordered_map<size_t, RegisteredFormat> index_to_parser;

//   template <typename T>
//   void register_format() {

//     size_t index = index_to_parser.size();
//     Metadata metadata = get_parser_metadata<T>();

//     index_to_parser.insert(
//         {index,
//          {.metadata = metadata,
//           .creator = [](std::string filepath, FileMode mode, FileCompression compression) {
//       return std::make_unique<T>(filepath, mode, compression);
//     }}});

//     // map extension to the index
//     for (const std::string& alias : metadata.aliases) {
//       std::unordered_map<std::string, size_t>::const_iterator it = alias_to_index.find(alias);
//       if (it != alias_to_index.end()) {
//         continue; // the alias is already used
//       }
//       alias_to_index.insert({alias, index});
//     }
//   }

//   RegisteredFormat from_alias(const std::string& alias) const {
//     std::unordered_map<std::string, size_t>::const_iterator it = alias_to_index.find(alias);
//     if (it == alias_to_index.end()) {
//       return {.metadata = Metadata{}, .creator = nullptr, .is_valid = false};
//     }
//     return index_to_parser.find(it->second)->second;
//   }

//   inline bool has_alias(const std::string& alias) const {
//     return alias_to_index.find(alias) != alias_to_index.end();
//   }

//   std::string infos() const;
// };

// const ParserFactory& parser_factory();

// /* ------------------------------------------------------------------------- */
// // Parser for LAMMPS Data format

// namespace lmpdata {

// enum class AtomStyle { ATOMIC, MOLECULAR, FULL, NONE };

// enum class Field { NATOM, TYPAT, XAXIS, YAXIS, ZAXIS, TILT, SECTION, END };

// enum class Section {
//   HEADER,
//   MASSES,
//   ATOMS,
//   VELOCITIES,
//   IGNORED,
// };

// constexpr std::array<Section, 3> sections_list { Section::MASSES, Section::ATOMS, Section::VELOCITIES };

// constexpr TokenNeedles<1> tok_atoms = {"atoms"};
// constexpr TokenNeedles<2> tok_atom_types = {"atom", "types"};
// constexpr TokenNeedles<2> tok_xaxis = {"xlo", "xhi"};
// constexpr TokenNeedles<2> tok_yaxis = {"ylo", "yhi"};
// constexpr TokenNeedles<2> tok_zaxis = {"zlo", "zhi"};
// constexpr TokenNeedles<3> tok_tilt = {"xy", "xz", "yz"};
// constexpr TokenNeedles<3> tok_sections = {"Masses", "Atoms", "Velocities"};
// constexpr TokenNeedles<3> tok_atom_style = {"atomic", "full", "molecule"};


// class LAMMPSDataParser : public TextParser {
// public:
//   LAMMPSDataParser(std::string filepath, FileMode mode, FileCompression compression)
//       : TextParser(filepath, mode, compression), m_atom_style(AtomStyle::NONE),
//         m_current_section(Section::HEADER) {}

//   int64_t next() override;
//   bool parse(Context& proxy) override;

// private:

//   AtomStyle m_atom_style;
//   Section m_current_section;
//   
//   // keep track of the current section line
//   inline std::string_view& section_line(void) { return line_buffer[1]; };
//   inline const std::string_view& section_line(void) const { return line_buffer[1]; };

//   bool parse_header(Context& ctx);
//   bool parse_atoms(Context& ctx);
//   bool parse_masses();
//   bool parse_velocities(Context& ctx);

//   void forward_to_next_section();
// };

// using AtomLineReader = bool(*)(const TokenSet&, Context&, size_t);
// using FieldParser = bool(*)(const TokenSet&, Context&);

// bool read_atom_line_atomic(const TokenSet&, Context&, size_t);

// inline const std::unordered_map<AtomStyle, AtomLineReader>& atom_line_readers() {
//   static std::unordered_map<AtomStyle, AtomLineReader> map = {
//     { AtomStyle::ATOMIC , &read_atom_line_atomic }
//   };
//   return map;
// }

// bool read_field_natom(const TokenSet&, Context&);
// bool read_field_atom_types(const TokenSet&, Context&);
// bool read_field_xaxis(const TokenSet&, Context&);
// bool read_field_yaxis(const TokenSet&, Context&);
// bool read_field_zaxis(const TokenSet&, Context&);
// bool read_field_tilt(const TokenSet&, Context&);
// bool read_field_section(const TokenSet&, Context&);

// inline const std::array<std::pair<Field, FieldParser>, 7>& field_readers() {
//   static constexpr std::array<std::pair<Field, FieldParser>, 7> readers = {{
//       {Field::NATOM, &read_field_natom},
//       {Field::TYPAT, &read_field_atom_types},
//       {Field::XAXIS, &read_field_xaxis},
//       {Field::YAXIS, &read_field_zaxis},
//       {Field::ZAXIS, &read_field_yaxis},
//       {Field::TILT, &read_field_tilt},
//       {Field::SECTION, &read_field_section},
//   }};
//   return readers;
// }

// inline Section get_section(std::string_view str) {
//   size_t index;
//   if (match_token_any(str, tok_sections, index))
//     return sections_list[index];
//   else
//     return Section::IGNORED;
// }

// } // namespace lmpdata

// template<> inline const Metadata get_parser_metadata<lmpdata::LAMMPSDataParser>() {
//   return {
//     .name        = "LAMMPS Data",
//     .aliases     = { "lmp", "lmp-data" },
//     .description = "LAMMPS Data format"
//   };
// }

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

// namespace xyz {

// enum class FieldType : uint8_t {
//   N = 0,
//   R = 1 << 0,
//   S = 1 << 1,
//   I = 1 << 2
// };

// inline constexpr FieldType operator|(FieldType lhs, FieldType rhs) {
//   return static_cast<FieldType>(static_cast<uint8_t>(lhs) | static_cast<uint8_t>(rhs));
// }

// inline constexpr FieldType operator&(FieldType lhs, FieldType rhs) {
//   return static_cast<FieldType>(static_cast<uint8_t>(lhs) & static_cast<uint8_t>(rhs));
// }

// inline constexpr FieldType operator~(FieldType f) { return static_cast<FieldType>(~static_cast<uint8_t>(f)); }

// inline constexpr FieldType& operator|=(FieldType& lhs, FieldType rhs) {
//   lhs = lhs | rhs;
//   return lhs;
// }

// inline constexpr FieldType& operator&=(FieldType& lhs, FieldType rhs) {
//   lhs = lhs & rhs;
//   return lhs;
// }

// enum class Field { SPECIES, POSITIONS, VELOCITIES, ID, TYPE, NIL };

// struct AtomicProperty {
//   Field field = Field::NIL;
//   FieldType type = FieldType::N;
//   size_t size = 0;
//   size_t index = 0;
// };

// struct Properties {
//   size_t index = 0;
//   std::array<AtomicProperty, 4> data{};

//   constexpr inline const AtomicProperty& species() const { return *(data.data() + 0); }
//   constexpr inline const AtomicProperty& positions() const { return *(data.data() + 1); }
//   constexpr inline const AtomicProperty& velo() const { return *(data.data() + 2); }
//   constexpr inline const AtomicProperty& id() const { return *(data.data() + 3); }

//   constexpr inline const bool has_species() const { return species().field == Field::SPECIES; }
//   constexpr inline const bool has_positions() const { return positions().field == Field::POSITIONS; }
//   constexpr inline const bool has_velo() const { return velo().field == Field::VELOCITIES; }
//   constexpr inline const bool has_id() const { return positions().field == Field::ID; }
// };

// // order matter - same order as in Properties
// constexpr TokenNeedles<4> tok_props{"species", "pos", "velo", "id"};
// constexpr std::array<Field, 4> props_list{Field::SPECIES, Field::POSITIONS, Field::VELOCITIES, Field::ID};

// constexpr TokenNeedles<3> tok_field_type = {"R", "S", "I"};
// constexpr std::array<FieldType, 3> type_list = {FieldType::R, FieldType::S, FieldType::I};

// template <size_t N1, size_t N2>
// bool read_exyz_comment_line(std::string_view line, TokenSetTmpl<N1>& tok_lattice, TokenSetTmpl<N2>& tok_props) {
//   constexpr std::string_view lattice_key = "Lattice=\"";
//   constexpr size_t lattice_key_len = 9;
//   constexpr std::string_view prop_key = "Properties=";
//   constexpr size_t prop_key_len = 11;

//   // constexpr IsSpace is_space{};
//   constexpr IsChar<'"'> is_quote{};
//   constexpr IsChar<' '> is_space{};

//   const char* begin = line.data();
//   const char* end = line.data() + line.size();

//   if (!starts_with(line, lattice_key))
//     return false;

//   const char* val_start = begin + lattice_key_len;
//   const char* val_end = val_start;

//   // while (val_end < end && *val_end != '"')
//   while (val_end < end && is_quote(*val_end))
//     ++val_end;
//   if (val_end >= end)
//     return false;

//   std::string_view lattice_view = std::string_view(val_start, val_end);
//   tokenize(lattice_view, tok_lattice);

//   const char* prop_it = val_end;
//   while (prop_it + prop_key_len <= end) {

//     if (std::memcmp(prop_it, prop_key.data(), prop_key_len) == 0) {
//       prop_it += prop_key_len;
//       const char* prop_end = prop_it;
//       // while (prop_end < end && *prop_end != ' ')
//       while (prop_end < end && is_space(prop_end))
//         ++prop_end;
//       std::string_view ppt = std::string_view(prop_it, prop_end);
//       tokenizer_tmpl(ppt, tok_props, ColumnDelimiter{});
//       break;
//     }

//     ++prop_it;
//   }

//   if (tok_props.size() == 0 || tok_props.size() % 3 != 0) return false;

//   return true;
// };

// inline bool parse_exyz_lattice(const TokenSet& tokens, IOContext& ctx) {
//   // read (ax ay az) (bx by bz) (cx cy cz)
//   return tokens_to_num(
//     tokens.at<0>(), ctx.H.m11,
//     tokens.at<1>(), ctx.H.m21,
//     tokens.at<2>(), ctx.H.m31,
//     tokens.at<3>(), ctx.H.m12,
//     tokens.at<4>(), ctx.H.m22,
//     tokens.at<5>(), ctx.H.m32,
//     tokens.at<6>(), ctx.H.m13,
//     tokens.at<7>(), ctx.H.m23,
//     tokens.at<8>(), ctx.H.m33
//   );
// }

// template <size_t N> inline bool parse_exyz_properties(const TokenSetTmpl<N>& tokens, Properties& properties) {
//   size_t begin = 0;
//   size_t size = 0;
//   for (size_t i = 0; i < tokens.size(); i += 3) {
//     size_t j = i + 1;
//     size_t k = i + 2;

//     tokens_to_num(tokens[k], size);

//     size_t index_field;
//     size_t index_type;

//     if (match_token_any(tokens[i], tok_props, index_field) && match_token_any(tokens[j], tok_field_type, index_type)) {
//       properties.data[index_field] = {
//           .field = props_list[index_field], .type = type_list[index_type], .size = size, .index = begin};
//     }
//     // TODO: validate type and size
//     begin += size;
//   }

//   return true;
// }

// inline bool parse_line(IOContext& ctx, const Properties& properties, const TokenSet& tokens, size_t particle_count) {

//   std::string_view type_as_string;
//   double x, y, z;
//   size_t type, id;

//   // parse type
//   if (!ctx.species.get_type_from_sv(tokens[properties.species().index], type)) {
//     lerr << "Undefined type: " << type_as_string << std::endl;
//     return false;
//   }

//   // parse positions
//   size_t index = properties.positions().index;
//   if (!(tokens_to_num(tokens[index], x, tokens[index + 1], y, tokens[index + 2], z))) {
//     lerr << "Unable to parse positions..." << std::endl;
//     return false;
//   }

//   // optional atom id
//   (properties.has_id()) ? tokens_to_num(tokens[properties.id().index], id) : id = particle_count;

//   // TODO: optional velocity

//   ParticleTupleIO& p = ctx.particle_data[particle_count];
//   p[field::rx] = x;
//   p[field::ry] = y;
//   p[field::rz] = z;
//   p[field::vx] = 0.0;
//   p[field::vy] = 0.0;
//   p[field::vz] = 0.0;
//   p[field::id] = id;
//   p[field::type] = type;

//   return true;
// }

// class ExtendedXYZParser : public TextParser {
// public:
//   ExtendedXYZParser(std::string filepath, FileMode mode, FileCompression compression)
//       : TextParser(filepath, mode, compression) {}

//   int64_t next() override;
//   bool parse(IOContext& ctx) override;

// private:
//   TokenSet64 m_ptokens;
//   Properties m_properties;
// };

// }

// template<> inline const Metadata get_parser_metadata<xyz::ExtendedXYZParser>() {
//   return {
//     .name        = "Extended XYZ",
//     .aliases     = { "xyz" },
//     .description = "Extended XYZ"
//   };
// }

// /* ------------------------------------------------------------------------- */

// struct FileInfos {
//   FileCompression compression = FileCompression::NONE;
//   const RegisteredFormat format{};
// };

// std::string remove_leading_chars(const std::string& str, const char c);
// const FileInfos guess_file_infos(const std::string& filepath);
// const FileInfos get_file_infos(const std::string& filepath, std::string format, std::string compression);

// } // namespace ioextra
// } // namespace exaStamp
