#pragma once

#include <array>
#include <cassert>
#include <cstdlib>
#include <cstring>

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

/* ------------------------------------------------------------------------- */
// Tokenizer utilities

using sv_span = std::span<const std::string_view>;

template <size_t Size>
struct TokenSetTmpl {

  using sv = std::string_view;

  alignas(64) std::array<sv, Size> tokens;
  size_t len = 0;

  // reset token count
  inline void clear() noexcept { len = 0; }

  inline sv* data() noexcept { return tokens.data(); }
  inline const sv* data() const noexcept { return tokens.data(); }

  inline sv& operator[](size_t i) noexcept {
    return tokens[i];
  }

  inline const sv& operator[](size_t i) const noexcept {
    return tokens[i];
  }

  template<size_t index>
  inline constexpr sv& at() noexcept {
    static_assert(index < Size);
    return tokens[index];
  }

  template<size_t index>
  inline constexpr const sv& at() const noexcept {
    static_assert(index < Size);
    return tokens[index];
  }

  constexpr size_t size() const noexcept { return len; }
  constexpr size_t capacity() const noexcept { return Size; }
  inline bool full() const noexcept { return len >= Size; }

  inline void push_back(sv token) noexcept {
    // TODO: Check if full ?
    tokens[len++] = token;
  } 

  // iteration accessor
  auto begin() const noexcept { return tokens.begin(); }
  auto end() const noexcept { return tokens.begin() + len; }
  auto begin() noexcept { return tokens.begin(); }
  auto end() noexcept { return tokens.begin() + len; }
};

using TokenSet = TokenSetTmpl<32>; // rather never exceeding 32 token per line
using TokenSet64 = TokenSetTmpl<64>; 

/* ------------------------------------------------------------------------- */

// Convert a strinv_view token to numeric type
template <typename T>
bool token_to_num(const std::string_view token, T& rval) {
  auto [p, rc] = std::from_chars(token.begin(), token.begin() + token.size(), rval);
  return rc == std::errc();
}

// Dark magic that expand to tokens_to_num with automatic type infering
#define parb ()
#define IOexpand(...) IOexpandB(IOexpandB(IOexpandB(IOexpandB(__VA_ARGS__))))
#define IOexpandB(...) IOexpandA(IOexpandA(IOexpandA(IOexpandA(__VA_ARGS__))))
#define IOexpandA(...) __VA_ARGS__

#define call_parse(token, value) token_to_num<decltype(value)>(token, value)
#define IOfor_each(macro, ...) __VA_OPT__(IOexpand(IOfor_each_helper(macro, __VA_ARGS__)))
#define IOfor_each_helper(macro, token, value, ...) macro(token, value) __VA_OPT__(IOfor_each_again parb (macro, __VA_ARGS__))
#define IOfor_each_again() && IOfor_each_helper
#define tokens_to_num(...) IOfor_each(call_parse, __VA_ARGS__)

/* ------------------------------------------------------------------------- */

inline bool is_space(char c) noexcept { return std::isspace(c); }

struct IsSpace {
  constexpr bool operator()(char c) const noexcept {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r';
  }
};

struct SpaceDelimiter {
  constexpr bool operator()(char c) const noexcept { return c == ' ' || c == '\t'; }
};

struct CommaDelimiter {
  constexpr bool operator()(char c) const noexcept { return c == ','; }
};

template <typename Predicate>
inline std::string_view trim(std::string_view str, Predicate&& is_trim_char) noexcept {
  while (!str.empty() && is_trim_char(str.front()))
    str.remove_prefix(1);
  while (!str.empty() && is_trim_char(str.back()))
    str.remove_suffix(1);
  return str;
}

template <typename Predicate>
inline std::string_view ltrim(std::string_view str, Predicate&& is_trim_char) noexcept {
  while (!str.empty() && is_trim_char(str.front()))
    str.remove_prefix(1);
  return str;
}

template <typename Predicate>
inline std::string_view rtrim(std::string_view str, Predicate&& is_trim_char) noexcept {
  while (!str.empty() && is_trim_char(str.back()))
    str.remove_suffix(1);
  return str;
}

inline std::string_view trim_spaces(std::string_view str) noexcept {
  return trim(str, IsSpace{});
}

// Old version of split function
// template <typename Callback, typename Predicated>
// inline void split(std::string_view str, Predicated&& is_delimiter, Callback&& callback) noexcept {
//   std::size_t pos = 0;
//   while (pos < str.size()) {
//     auto const next_pos = std::find_if(str.begin() + pos, str.end(), is_delimiter) - str.begin();
//     callback(str.substr(pos, next_pos - pos));
//     pos = static_cast<std::size_t>(next_pos) == str.size() ? str.size() : next_pos + 1;
//   }
// }

// new split function
template <typename Callback, typename Delimiter>
inline void split(std::string_view str, Delimiter&& is_delimiter, Callback&& callback) noexcept {
  const char* data = str.data();
  const char* end = data + str.size();
  while (data < end) {
    while (data < end && is_delimiter(*data))
      ++data;
    const char* token_start = data;
    while (data < end && !is_delimiter(*data))
      ++data;
    if (token_start != data)
      callback(std::string_view(token_start, data - token_start));
  }
}

//TODO: tokenize function that tokenize up to a maximum number of token

template <typename Delimiter>
inline void tokenizer_tmpl(const std::string_view str, TokenSet& tokens, Delimiter&& delimiter) {
  tokens.clear(); // reset the token set
  split(str, delimiter, [&tokens](std::string_view token) {
    if (!token.empty() && token.front() != '#') {
      tokens.push_back(token);
    }
  });
}

inline void tokenize(const std::string_view str, TokenSet& tokens) {
  return tokenizer_tmpl(str, tokens, SpaceDelimiter{});
}


template<size_t N> using TokenNeedles = std::array<std::string_view, N>;

// return true at the first needle that the given token
inline bool _token_match_any(std::string_view t, sv_span n) {
  const size_t nc = n.size();
  if (nc == 0)
    return false;
  for (size_t i = 0; i < nc; ++i)
    if (t == n[i])
      return true;
  return false;
}

// return true at the first needle that the given token with the needle index
inline bool _token_match_any(std::string_view t, sv_span n, size_t& in) {
  const size_t nc = n.size();
  if (nc == 0)
    return false;
  for (size_t i = 0; i < nc; ++i)
    if (t == n[i]) {
      in = i;
      return true;
    }
  return false;
}

// return true at the first needle that match a token in the set
inline bool _token_set_match_any(sv_span t /*tokens*/, sv_span n /*needles*/) {
  const size_t tc = t.size();
  const size_t nc = n.size();

  if (nc == 0)
    return false;

  for (size_t i = 0; i < tc; ++i)
    for (size_t j = 0; j < nc; ++j)
      if (t[i] == n[j])
        return true;
  return false;
}

// return true at the first needle that match a token in the set
inline bool _token_set_match_any(sv_span t /*tokens*/, sv_span n /*needles*/, size_t& it, size_t in) {
  const size_t tc = t.size();
  const size_t nc = n.size();

  if (nc == 0)
    return false;

  for (size_t i = 0; i < tc; ++i)
    for (size_t j = 0; j < nc; ++j)
      if (t[i] == n[j]) {
        it = i; // token index
        in = j; // needle index
        return true;
      }
  return false;
}

inline bool _token_match_sequence(sv_span t, sv_span n) {
  const size_t tc = t.size();
  const size_t nc = n.size();

  if (nc == 0 || tc < nc) return false;

  for (size_t i = 0; i <= tc - nc; ++i) {
      bool match = true;
      for (size_t j = 0; j < nc; ++j) {
          if (t[i + j] != n[j]) {
              match = false;
              break;
          }
      }
      if (match) return true;
  }
  return false;
}

template <size_t N>
inline bool match_token_any(const std::string_view token, const TokenNeedles<N>& needles) {
  return _token_match_any(token, sv_span(needles.data(), N));
}

template <size_t N>
inline bool match_token_any(const std::string_view token, const TokenNeedles<N>& needles, size_t& i) {
  return _token_match_any(token, sv_span(needles.data(), N), i);
}

template <size_t N>
inline bool match_token_set_any(const TokenSet& tokens, const TokenNeedles<N>& needles) {
  return _token_set_match_any(sv_span(tokens.data(), tokens.size()), sv_span(needles.data(), N));
}

template <size_t N>
inline bool match_token_set_any(const TokenSet& tokens, const TokenNeedles<N>& needles, size_t& i, size_t& j) {
  return _token_set_match_any(sv_span(tokens.data(), tokens.size()), sv_span(needles.data(), N), i, j);
}

template <size_t N>
inline bool match_token_set_any(const TokenSet& tokens, const TokenNeedles<N>& needles, size_t& i) {
  size_t j = 0;
  return _token_set_match_any(sv_span(tokens.data(), tokens.size()), sv_span(needles.data(), N), i, j);
}

template <size_t N>
inline bool match_token_sq(const TokenSet& tokens, const TokenNeedles<N>& needles) {
  return _token_match_sequence(sv_span(tokens.data(), tokens.size()), sv_span(needles.data(), N));
}

/* ------------------------------------------------------------------------- */

// namespace re {

// using match_str = std::match_results<std::string::const_iterator>;
// using match_strv = std::match_results<std::string_view::const_iterator>;

// // regex iterator for string view
// using sv_regex_iterator = std::regex_iterator<std::string_view::const_iterator>;

// template <typename T> struct RegexMatches {

//   inline T& operator[](size_t i) {
//     assert(i < m_data.size());
//     return m_data[i];
//   }

// private:
//   std::array<T, 10> m_data;
// };

// using RegexMatchesStringView = RegexMatches<match_strv>;

// inline size_t regex_find_all(const char* begin, const char* end, const std::regex& re_pattern) {

//   size_t count = 0;
//   sv_regex_iterator it_begin = sv_regex_iterator(begin, end, re_pattern);
//   sv_regex_iterator it_end = sv_regex_iterator();

//   for (; it_begin != it_end; ++it_begin) {
//     ++count;
//   }
//   return count;
// };

// static const std::string str_word = std::string("([a-zA-Z]+)");
// static const std::string str_uint = std::string("([0-9]+)");
// static const std::string str_int = std::string("([+-]?[0-9]+)");
// static const std::string str_float = std::string("([+-]?\\d*[.]?\\d*[Eefd]?[+-]?\\d+)");
// static const std::string str_any = std::string("(\\S+)");
// static const std::string str_space = std::string("\\s+");
// static const std::string str_lstart_f = std::string("^(?!#)\\s*") + str_float;
// static const std::string str_lstart_i = std::string("^(?!#)\\s*") + str_int;
// static const std::string str_vec2d = str_float + str_space + str_float;
// static const std::string str_vec3d = str_float + str_space + str_float + str_space + str_float;
// static const std::string str_mat3d = str_vec3d + str_space + str_vec3d + str_space + str_vec3d;

// const std::regex re_uint = std::regex(str_uint);
// const std::regex re_float = std::regex(str_float);
// const std::regex re_lsart_f = std::regex(str_lstart_f);
// const std::regex re_lstart_i = std::regex(str_lstart_i);
// const std::regex re_vec2d = std::regex(str_vec2d);
// const std::regex re_vec3d = std::regex(str_vec3d);
// const std::regex re_mat3d = std::regex(str_mat3d);

// } // namespace re

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
  Mat3d M = Mat3d{};                // transpose(H)
  Mat3d invM = Mat3d{};             // Ht^-1
  Mat3d Xform = Mat3d{};            // initial xform (likely == I)
  Mat3d invXform = Mat3d{};         // intial xform^-1
  Mat3d domainXform = Mat3d{};      // new xform
  Vec3d boxlen{};                   // box len in x, y, z
  bool uniform_scale;               // ?

  // for timing purpose
  size_t elapsed_time_ms = 0;
  size_t elapsed_time_ys = 0;
  size_t elapsed_time_ns = 0;

  // TODO: Move this out
  inline void init_domain(Domain& domain, ReadBoundsSelectionMode bounds_mode) {

    Vec3d a{H.m11, H.m12, H.m13};
    Vec3d b{H.m21, H.m22, H.m23};
    Vec3d c{H.m31, H.m32, H.m33};

    boxlen.x = norm(a);
    boxlen.y = norm(b);
    boxlen.z = norm(c);

    lout << "a = " << a << onika::format_string(" , norm = %10.5f ", boxlen.x) << std::endl;
    lout << "b = " << b << onika::format_string(" , norm = %10.5f ", boxlen.y) << std::endl;
    lout << "c = " << c << onika::format_string(" , norm = %10.5f ", boxlen.z) << std::endl;

    if (!domain.xform_is_identity()) {
      lout << std::endl << "Init initial xform to identity" << std::endl;
      domain.set_xform(make_identity_matrix());
    }

    AABB domain_bounds = {{0., 0., 0.}, boxlen};
    compute_domain_bounds(domain, bounds_mode, 0.0, domain_bounds, domain_bounds, true);

    M = transpose(H);
    invM = inverse(M);

    invXform = (domain.xform_is_identity()) ? make_identity_matrix() : domain.inv_xform();
    domainXform = M * inverse(diag_matrix(boxlen)) * domain.xform();
    uniform_scale = is_uniform_scale(domainXform);

    if (uniform_scale) {
      domain.set_xform(make_identity_matrix());
      domain.set_cell_size(domain.cell_size() * domainXform.m11);
      domain.set_bounds({domain.origin() * domainXform.m11, domain.extent() * domainXform.m11});
    } else {
      domain.set_xform(domainXform);
    }
  }
};

using IOContext = Context;

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
  TokenSet m_tokens;

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

// namespace xyz {
// using namespace re;

// static const std::string str_lattice = std::string("Lattice=.") + str_mat3d; // we accept any quote "' or none
// static const std::string str_props = std::string("Properties=\\S+");
// static const std::string str_prop = std::string("([a-z]+:[SR]:[0-9])");
// static const std::string str_pp = std::string("([^:]+):([^:]+):([^:]+)");

// static const std::regex re_lattice = std::regex(str_lattice, std::regex_constants::icase);
// static const std::regex re_props = std::regex(str_props, std::regex_constants::icase);
// static const std::regex re_prop = std::regex(str_prop);
// static const std::regex re_pp = std::regex(str_pp);

// enum class Type : uint8_t {
//   N = 0,
//   R = 1 << 0,
//   S = 1 << 1,
//   I = 1 << 2
// };

// inline constexpr Type operator|(Type lhs, Type rhs) {
//   return static_cast<Type>(static_cast<uint8_t>(lhs) | static_cast<uint8_t>(rhs));
// }

// inline constexpr Type operator&(Type lhs, Type rhs) {
//   return static_cast<Type>(static_cast<uint8_t>(lhs) & static_cast<uint8_t>(rhs));
// }

// inline constexpr Type operator~(Type f) { return static_cast<Type>(~static_cast<uint8_t>(f)); }

// inline constexpr Type& operator|=(Type& lhs, Type rhs) {
//   lhs = lhs | rhs;
//   return lhs;
// }

// inline constexpr Type& operator&=(Type& lhs, Type rhs) {
//   lhs = lhs & rhs;
//   return lhs;
// }

// enum class Field {
//   SPECIES,
//   POSITIONS,
//   VELOCITIES,
//   ID,
//   TYPE,
//   NIL
// };

// using FieldMap = std::unordered_map<std::string, Field>;
// using TypeMap = std::unordered_map<std::string, Type>;
// using ValidatorMap = std::unordered_map<Field, std::pair<Type, size_t>>;

// inline const FieldMap& field_map() {
//   static const FieldMap map = {
//     {"species", Field::SPECIES},
//     {"pos"    , Field::POSITIONS},
//     {"id"     , Field::ID},
//     {"vel"    , Field::VELOCITIES},
//     {"type"   , Field::TYPE},
//     {"nil"    , Field::NIL}
//   };
//   return map;
// }

// inline const TypeMap& type_map() {
//   static const TypeMap map = {
//       {"R", Type::R},
//       {"S", Type::S},
//       {"I", Type::I},
//   };
//   return map;
// }

// inline const std::unordered_map<Type, std::string>& type_regex() {
//   static const std::unordered_map<Type, std::string> map = {
//       {Type::R, str_float},
//       {Type::S, str_word},
//       {Type::I, str_int},
//   };
//   return map;
// }

// inline const ValidatorMap& validators() {
//   static const ValidatorMap map = {
//     {Field::SPECIES   , {Type::S, 1}},
//     {Field::POSITIONS , {Type::R, 3}},
//     {Field::ID        , {Type::I, 1}},
//     {Field::VELOCITIES, {Type::R, 3}},
//     {Field::TYPE      , {Type::R | Type::I, 1} },
//   };
//   return map;
// }

// struct Property {
//   Field field = Field::NIL;
//   Type type = Type::N;
//   size_t size = 0;
//   size_t begin = 0;
// };

// struct Properties {
//   size_t index = 0;
//   std::vector<Property> properties;
//   std::unordered_map<Field, size_t> indices;

//   bool add(const std::string&, const std::string&, size_t);
//   std::string build_line_regex_pattern(void);
//   inline Property& get(Field f) { return properties[indices.find(f)->second]; };
// };

// class ExtendedXYZParser : public TextParser {
// public:
//   ExtendedXYZParser(std::string filepath, FileMode mode, FileCompression compression)
//       : TextParser(filepath, mode, compression) {}

//   int64_t next() override;
//   bool parse(Context& ctx) override;
// };

// } // namespace xyz

// template<> inline const Metadata get_parser_metadata<xyz::ExtendedXYZParser>() {
//   return {
//     .name        = "Extended XYZ",
//     .aliases     = { "xyz" },
//     .description = "Extended XYZ"
//   };
// }

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
