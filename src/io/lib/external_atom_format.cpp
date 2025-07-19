#include <filesystem>
#include <exaStamp/io/external_atom_format.h>

namespace exaStamp {
using namespace exanb;

namespace ioextra {

/* ------------------------------------------------------------------------- */

TextFileHandler::TextFileHandler(const std::string& filepath) : m_path(filepath) {}

ASCIIFileHandler::ASCIIFileHandler(const std::string& path, FileMode mode) : TextFileHandler(path) {
  m_fptr = std::fopen(m_path.c_str(), convert_mode_to_char(mode));
  if (m_fptr == nullptr) {
    std::abort();
  }
}

ASCIIFileHandler::~ASCIIFileHandler() {
  if (m_fptr != nullptr)
    std::fclose(m_fptr);
}

void ASCIIFileHandler::clear() noexcept { std::clearerr(m_fptr); }

void ASCIIFileHandler::seek(uint64_t cursor) {
  static_assert(sizeof(uint64_t) == sizeof(off64_t));

  auto status = fseeko(m_fptr, static_cast<off64_t>(cursor), SEEK_SET);
  if (status != 0) {
    // auto message = std::strerror(errno);
    std::abort();
  }
}

size_t ASCIIFileHandler::read(char* data, size_t size) {
  size_t result = std::fread(data, 1, size, m_fptr);
  if (std::ferror(m_fptr) != 0)
    std::abort();
  return result;
}

#ifdef USE_ZLIB

GzipFileHandler::GzipFileHandler(const std::string& path, FileMode mode) : TextFileHandler(path) {
  const char* openmode;
  switch (mode) {
  case FileMode::READ:
    openmode = "rb";
    break;
  case FileMode::APPEND:
  case FileMode::WRITE:
    std::abort();
    break;
  }

  m_fptr = gzopen64(m_path.c_str(), openmode);
  if (m_fptr == nullptr) {
    lerr() << onika::format_string("Unable to open file: ", m_path.c_str()) << std::endl;
    std::abort();
  }
}

GzipFileHandler::~GzipFileHandler() {
  if (m_fptr != nullptr)
    gzclose(m_fptr);
}

void GzipFileHandler::clear() noexcept { gzclearerr(m_fptr); }

void GzipFileHandler::seek(uint64_t cursor) {
  static_assert(sizeof(uint64_t) == sizeof(z_off64_t));

  auto status = gzseek64(m_fptr, static_cast<z_off64_t>(cursor), SEEK_SET);
  if (status == -1) {
    const char* message = gz_error();
    lerr() << onika::format_string("Error while decompressing gzip file : %s", message) << std::endl;
    std::abort();
  }
}

size_t GzipFileHandler::read(char* buffer, size_t size) {
  int read = gzread(m_fptr, buffer, safe_cast(size));
  const char* error = gz_error();
  if (read == -1 || error != nullptr) {
    lerr() << onika::format_string("Error while reading gzip file : %s", error) << std::endl;
    std::abort();
  }
  return static_cast<size_t>(read);
}

const char* GzipFileHandler::gz_error(void) const {
  int status = Z_OK;
  const char* message = gzerror(m_fptr, &status);
  return status == Z_OK ? nullptr : message;
}

unsigned GzipFileHandler::safe_cast(size_t value) {
  constexpr size_t max = std::numeric_limits<unsigned>::max();
  if (value > max) {
    lerr() << onika::format_string("%d is too big for unsigned in call to zlib function", value) << std::endl;
    std::abort();
  }
  return static_cast<unsigned>(value);
}

#endif

#ifdef USE_LZMA

size_t XzFileHandler::safe_cast(size_t size) {
  constexpr size_t max = std::numeric_limits<size_t>::max();
  if (size > max) {
    lerr() << onika::format_string("%ld is too big for unsigned in call to zlib function", size) << std::endl;
    std::abort();
  }
  return static_cast<size_t>(size);
}

void XzFileHandler::check_lzma_ret(lzma_ret ret) {
  int retnum = static_cast<int>(ret);

  switch (ret) {
  case LZMA_OK:
  case LZMA_STREAM_END:
    return;
  case LZMA_GET_CHECK:
  case LZMA_NO_CHECK:
    lerr() << onika::format_string("lzma: no integrity check.", retnum) << std::endl;
    std::abort();
  case LZMA_MEM_ERROR:
  case LZMA_MEMLIMIT_ERROR:
    lerr() << onika::format_string("lzma: failed to allocated memory (err: %d)", retnum) << std::endl;
    std::abort();
  case LZMA_FORMAT_ERROR:
    lerr() << onika::format_string("lzma: invalid .xz format (err: %d)", retnum) << std::endl;
    std::abort();
  case LZMA_OPTIONS_ERROR:
    lerr() << onika::format_string("lzma: invalid options (err: %d)", retnum) << std::endl;
    std::abort();
  case LZMA_DATA_ERROR:
  case LZMA_BUF_ERROR:
    lerr() << onika::format_string("lzma: file is corrupted or truncated (err: %d)", retnum) << std::endl;
    std::abort();
  case LZMA_UNSUPPORTED_CHECK:
    lerr() << onika::format_string("lzma: file intergrity check not supported (err: %d)", retnum) << std::endl;
    std::abort();
  case LZMA_PROG_ERROR:
    lerr() << onika::format_string("lzma: this is bug (err: %d)", retnum) << std::endl;
    std::abort();
  }
}

void XzFileHandler::start_lzma_decoder_stream(lzma_stream* stream) {
  constexpr uint64_t memory_limit = std::numeric_limits<uint64_t>::max();
  check_lzma_ret(lzma_stream_decoder(stream, memory_limit, LZMA_TELL_UNSUPPORTED_CHECK | LZMA_CONCATENATED));
}

XzFileHandler::XzFileHandler(const std::string& path, FileMode mode)
    : TextFileHandler(path), m_xz_buffer(IOEXTRA_LZMA_INTERNAL_BUF_SIZE) {
  const char* openmode;
  switch (mode) {

  case FileMode::READ:
    openmode = "rb";
    start_lzma_decoder_stream(&m_xz_stream);
    break;

  case FileMode::APPEND:
  case FileMode::WRITE:
    lerr() << "Only reading is suported." << std::endl;
    std::abort();
    break;
  }

  m_fptr = std::fopen(m_path.c_str(), openmode);
  if (m_fptr == nullptr) {
    lzma_end(&m_xz_stream);
    lerr() << onika::format_string("Unable to open file: ", m_path.c_str()) << std::endl;
    std::abort();
  }
}

XzFileHandler::~XzFileHandler() {
  lzma_end(&m_xz_stream);
  if (m_fptr != nullptr)
    std::fclose(m_fptr);
}

void XzFileHandler::clear() noexcept { std::clearerr(m_fptr); }

void XzFileHandler::seek(uint64_t cursor) {
  // lzma is a stream based compression format, so random access to file is not supported.
  // inefficient implementation: re-decompressing the file from the begining
  // not really an issue if we read from the begining of the file.

  lzma_end(&m_xz_stream);
  m_xz_stream = LZMA_STREAM_INIT;
  start_lzma_decoder_stream(&m_xz_stream);

  // set position to 0
  std::fseek(m_fptr, 0, SEEK_SET);
  constexpr size_t bufsize = IOEXTRA_LZMA_STREAM_BUF_SIZE;
  char buffer[bufsize];

  while (cursor > bufsize) {
    size_t read_size = this->read(buffer, bufsize);
    assert(read_size == bufsize);
    cursor -= read_size;
  }

  [[maybe_unused]] size_t read_size = this->read(buffer, static_cast<size_t>(cursor));
  assert(read_size == cursor);
}

size_t XzFileHandler::read(char* buffer, size_t size) {

  m_xz_stream.next_out = reinterpret_cast<uint8_t*>(buffer);
  m_xz_stream.avail_out = size;

  while (m_xz_stream.avail_out > 0) {
    if (m_xz_stream.avail_in == 0) {

      m_xz_stream.next_in = m_xz_buffer.data();
      m_xz_stream.avail_in = std::fread(m_xz_buffer.data(), 1, m_xz_buffer.size(), m_fptr);

      if (std::ferror(m_fptr)) {
        lerr() << "Something went wrong while reading xz file" << std::endl;
        std::abort();
      }
    }

    lzma_action action = std::feof(m_fptr) && m_xz_stream.avail_in == 0 ? LZMA_FINISH : LZMA_RUN;
    lzma_ret ret = lzma_code(&m_xz_stream, action);

    if (ret == LZMA_STREAM_END) {
      break;
    }

    check_lzma_ret(ret);
  }

  return size - m_xz_stream.avail_out;
}

#endif

#ifdef USE_BZIP2

unsigned Bzip2FileHandler::safe_cast(uint64_t size) {
  constexpr size_t max = std::numeric_limits<size_t>::max();
  if (size > max) {
    lerr() << onika::format_string("%ld is too big for unsigned in call to bzlib function", size) << std::endl;
    std::abort();
  }
  return static_cast<unsigned>(size);
}

void Bzip2FileHandler::check_bz2_retcode(int code) {

  switch (code) {
  case BZ_OK:
  case BZ_RUN_OK:
  case BZ_FLUSH_OK:
  case BZ_FINISH_OK:
  case BZ_STREAM_END:
    return;
  case BZ_SEQUENCE_ERROR:
  case BZ_PARAM_ERROR:
    lerr() << onika::format_string("bzip2: bad call to bzlib (err: %d)", code) << std::endl;
    std::abort();
  case BZ_MEM_ERROR:
    lerr() << onika::format_string("bzip2: memory allocation failed (err: %d)", code);
    std::abort();
  case BZ_DATA_ERROR:
    lerr() << onika::format_string("bzip2: corrupted file (err: %d)", code);
    std::abort();
  case BZ_DATA_ERROR_MAGIC:
    lerr() << onika::format_string("bzip2: this file do not seems to be a bz2 file (err: %d)", code);
    std::abort();
  // These errors should not occur when using the stream API
  case BZ_CONFIG_ERROR:
    lerr() << onika::format_string("bzip2: mis-compiled bzlib (err: %d)", code);
    std::abort();
  case BZ_IO_ERROR:
  case BZ_UNEXPECTED_EOF:
  case BZ_OUTBUFF_FULL:
    lerr() << onika::format_string("bzip2: unexpected error from bzlib (err: %d)", code);
    std::abort();
  default:
    lerr() << onika::format_string("bzip2: ???? (err: %d)", code);
    std::abort();
  }
}

Bzip2FileHandler::Bzip2FileHandler(const std::string& path, FileMode mode)
    : TextFileHandler(path), m_bz2_buffer(IOEXTRA_BZ2_INTERNAL_BUF_SIZE) {

  const char* openmode = nullptr;
  switch (mode) {
  case FileMode::READ:
    openmode = "rb";
    m_end_bz2_stream = BZ2_bzDecompressEnd;
    std::memset(&m_bz2_stream, 0, sizeof(bz_stream));
    check_bz2_retcode(BZ2_bzDecompressInit(&m_bz2_stream, 0, 0));
    break;

  case FileMode::APPEND:
  case FileMode::WRITE:
    lerr() << "Only reading is suported with bz2 compression format" << std::endl;
    std::abort();
    break;
  }

  m_fptr = std::fopen(m_path.c_str(), openmode);
  if (m_fptr == nullptr) {
    m_end_bz2_stream(&m_bz2_stream);
    lerr() << onika::format_string("Unable to open file: ", m_path.c_str()) << std::endl;
    std::abort();
  }
}

Bzip2FileHandler::~Bzip2FileHandler() {
  m_end_bz2_stream(&m_bz2_stream);
  if (m_fptr != nullptr)
    std::fclose(m_fptr);
}


void Bzip2FileHandler::clear() noexcept { std::clearerr(m_fptr); }

void Bzip2FileHandler::seek(uint64_t cursor) {
  // bzip2 is a stream based compression format, so random access to file is not supported.
  // inefficient implementation: re-decompressing the file from the begining
  // not really an issue if we read from the begining of the file.

  m_end_bz2_stream(&m_bz2_stream);
  std::memset(&m_bz2_stream, 0, sizeof(bz_stream));
  check_bz2_retcode(BZ2_bzDecompressInit(&m_bz2_stream, 0, 0));

  // set position to 0
  std::fseek(m_fptr, 0, SEEK_SET);
  constexpr size_t bufsize = IOEXTRA_BZ2_STREAM_BUF_SIZE;
  char buffer[bufsize];

  while (cursor > bufsize) {
    size_t read_size = read(buffer, bufsize);
    assert(read_size == bufsize);
    cursor -= read_size;
  }

  [[maybe_unused]] size_t read_size = read(buffer, static_cast<size_t>(cursor));
  assert(read_size == cursor);
}

size_t Bzip2FileHandler::read(char* buffer, size_t size) {

  m_bz2_stream.next_out = buffer;
  m_bz2_stream.avail_out = safe_cast(size);

  while (m_bz2_stream.avail_out > 0) {
    if (m_bz2_stream.avail_in == 0) {

      m_bz2_stream.next_in = m_bz2_buffer.data();
      m_bz2_stream.avail_in = safe_cast(std::fread(m_bz2_buffer.data(), 1, m_bz2_buffer.size(), m_fptr));

      if (std::ferror(m_fptr)) {
        lerr() << "Something went wrong while reading xz file" << std::endl;
        std::abort();
      }
    }

    int retcode = BZ2_bzDecompress(&m_bz2_stream);

    if (retcode == BZ_STREAM_END) {
      break;
    }
    // check for error
    check_bz2_retcode(retcode);
  }

  return size - m_bz2_stream.avail_out;
}

#endif

/* ------------------------------------------------------------------------- */

TextFile::TextFile(std::string filepath, FileMode mode, FileCompression compression)
    : File(std::move(filepath), mode, compression), m_buf(IOEXTRA_FILE_BUFFER_SIZE, 0), m_lstart(m_buf.data()),
      m_bend(m_buf.data() + m_buf.size()) {

  switch (compression) {
  case FileCompression::BZIP2:
#ifdef USE_BZIP2
    m_handler = std::make_unique<Bzip2FileHandler>(m_path, m_mode);
#else
    lerr() << "BZ2 Compression format is not supported" << std::endl;
    lerr() << "Compile with BZ2 support -DWITH_BZIP2=ON" << std::endl;
#endif
    break;

  case FileCompression::XZ:
#ifdef USE_LZMA
    m_handler = std::make_unique<XzFileHandler>(m_path, m_mode);
#else
    lerr() << "XZ Compression format is not supported" << std::endl;
    lerr() << "Compile with XZ support -DWITH_LZMA=ON" << std::endl;
#endif
    break;

  case FileCompression::GZIP:
#ifdef USE_ZLIB
    m_handler = std::make_unique<GzipFileHandler>(m_path, m_mode);
#else
    lerr() << "GZ Compression format is not supported" << std::endl;
    lerr() << "Compile with GZ support -DWITH_ZLIB=ON" << std::endl;
#endif
    break;

  case FileCompression::NONE:
    m_handler = std::make_unique<ASCIIFileHandler>(m_path, m_mode);
    break;
  }

  if (m_handler == nullptr) {
    lerr() << "Unable to instantiate a file handler." << std::endl;
    std::abort();
  }
}

uint64_t TextFile::tell() const {
  // assert that we are inside the buffer
  assert(m_lstart >= m_buf.data());
  uint64_t delta = is_buffer_init() ? static_cast<uint64_t>(m_lstart - m_buf.data()) : 0;
  return m_cursor + delta;
}

void TextFile::seek(uint64_t pos) {
  m_eof = false;
  m_hdlr_eof = false;

  if (is_buffer_init()) {
    int64_t delta = static_cast<int64_t>(pos) - static_cast<int64_t>(m_cursor);
    if (delta >= 0 && delta < static_cast<int64_t>(m_buf.size())) {
      m_lstart = m_buf.data() + delta;
      m_eof = false;
      return;
    }
  }

  m_handler->seek(pos);
  m_cursor = pos;
  m_buf[0] = IOEXTRA_NULL_CHAR;
}

void TextFile::clear() {
  m_eof = false;
  m_hdlr_eof = false;
  m_handler->clear();
}

void TextFile::reset() {
  clear();
  seek(0);
}

std::string_view TextFile::get_line() {

  // if the buffer is not initialized fill it at pos = 0
  if (!is_buffer_init())
    fill_buffer(0);

  if (m_eof)
    return "";

  size_t length = 0;
  size_t msvc = 0; // windows compatibility

  while (true) {

    // number of charater still in the buffer
    size_t remain = static_cast<size_t>(m_bend - m_lstart);

    // needle is the position of the next new line
    // const char* newline = reinterpret_cast<const char*>(std::memchr(m_lstart + length, '\n', remain - length));
    const char* newline = static_cast<const char*>(std::memchr(m_lstart + length, '\n', remain - length));

    if (newline != nullptr) {
      assert(m_lstart <= newline);
      length += static_cast<size_t>(newline - m_lstart + 1);

      // check for windows string compatibility
      if (newline > m_buf.data() && newline[-1] == '\r')
        msvc = 1;

      break;

    } else if (m_hdlr_eof) {
      // we reach the end of the file
      m_eof = true;

      if (m_lstart != m_bend - 1) {
        std::string_view line = std::string_view(m_lstart);
        m_lstart += line.length();
        return line;
      }
    }

    if (remain >= m_buf.size()) {
      size_t shift = m_lstart - m_buf.data();
      m_buf.resize(2 * m_buf.size(), 0);
      m_lstart = m_buf.data() + shift;
      m_bend = m_buf.data() + m_buf.size();
    }

    std::memmove(m_buf.data(), m_lstart, remain);
    fill_buffer(remain);
  }

  std::string_view line = std::string_view(m_lstart, length - msvc - 1);
  m_lstart += length;

  return line;
}

void TextFile::fill_buffer(size_t pos) {
  size_t size = m_buf.size() - pos;

  if (is_buffer_init())
    m_cursor += size;

  size_t read_size = m_handler->read(m_buf.data() + pos, size);

  if (read_size < size) {
    m_hdlr_eof = true;
    std::memset(m_buf.data() + pos + read_size, 0, size - read_size);
  }

  m_lstart = m_buf.data();
}

/* ------------------------------------------------------------------------- */

void TextParser::scan() {
  // if end of file, return.
  if (m_eof)
    return;

  size_t i = m_file.tell();

  while (!m_file.eof()) {
    lout << "kkkk" << std::endl;
    if (int64_t c = next(); c > -1) {
      loc.push_back(c);
    }
    break;
  }

  m_eof = true;
  m_file.clear();

  if (i == 0 && !loc.empty()) {
    m_file.seek(loc[0]);
  } else {
    m_file.seek(i);
  }
}

bool TextParser::operator()(Context& ctx) {
  std::chrono::time_point time_start = std::chrono::high_resolution_clock::now();

  m_file.seek(0);
  bool parsed = parse(ctx);

  std::chrono::time_point time_end = std::chrono::high_resolution_clock::now();
  ctx.elapsed_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
  ctx.elapsed_time_ys = std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
  ctx.elapsed_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_start).count();
  return parsed;
}

bool TextParser::operator()(Context& ctx, size_t index) {

  std::chrono::time_point time_start = std::chrono::high_resolution_clock::now();

  if (index >= size())
    scan();

  if (size() == 0) {
    lerr() << "File does not contain any valid step to read.\n" << std::endl;
    return false;
  }

  if (index >= size()) {
    lerr() << onika::format_string("Cant read step=%ld. Maximal step=%ld", index, size()) << std::endl;
    return false;
  }

  m_file.seek(static_cast<uint64_t>(loc[index]));
  bool parsed = parse(ctx);

  std::chrono::time_point time_end = std::chrono::high_resolution_clock::now();
  ctx.elapsed_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
  ctx.elapsed_time_ys = std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_start).count();
  ctx.elapsed_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_start).count();
  return parsed;
}

/* ------------------------------------------------------------------------- */

ParserFactory::ParserFactory() {

  register_format<lmpdata::LAMMPSDataParser>();
  // register_format<LAMMPSDumpParser>();
  register_format<xyz::ExtendedXYZParser>();
};

std::string ParserFactory::infos() const {

  std::string infos = "Supported atom format:\n";

  std::unordered_map<size_t, RegisteredFormat>::const_iterator it;
  size_t size = index_to_parser.size();

  for (size_t index = 0; index < size; ++index) {
    infos += onika::format_string("\t- %s\n", index_to_parser.find(index)->second.metadata.name.c_str());
  }

  infos += "\n";

  return infos;
};

const ParserFactory& parser_factory() {
  static const ParserFactory parser_factory_instance;
  return parser_factory_instance;
};

/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Data format

namespace lmpdata {

// LAMMPS Data only contains one frame
int64_t LAMMPSDataParser::next() {
  if (size_t cursor = m_file.tell(); cursor == 0) {
    m_file.get_line();
    return static_cast<int64_t>(cursor);
  } else {
    return -1;
  }
}


bool LAMMPSDataParser::parse(Context& ctx) {

  lout << "Parsing LAMMPS Data file (v2):" << std::endl;

  current_line() = m_file.get_line();

  while (!m_file.eof()) {

    current_line() = trim_spaces(current_line());
    if (current_line().empty() || current_line().front() == '#') {
      current_line() = m_file.get_line();
      continue;
    }

    switch (m_current_section) {
    case Section::HEADER:
      parse_header(ctx);
      break;
    case Section::MASSES:
      parse_masses();
      break;
    case Section::ATOMS:
      parse_atoms(ctx);
      break;
    case Section::IGNORED:
      forward_to_next_section();
      break;
    default:
      break;
    }

    current_line() = m_file.get_line();
  }

  // print some info after parsing
  lout << onika::format_string("- %-15s = %ld", "particle_count", ctx.n_particles) << std::endl;
  lout << onika::format_string("- %-15s = %ld", "atom_types", ctx.n_species) << std::endl;
  lout << onika::format_string("- %-15s = (%.5e, %.5e)", "(xlo, xhi)", ctx.origin.x, ctx.origin.x + ctx.H.m11) << std::endl;
  lout << onika::format_string("- %-15s = (%.5e, %.5e)", "(ylo, yhi)", ctx.origin.y, ctx.origin.y + ctx.H.m22) << std::endl;
  lout << onika::format_string("- %-15s = (%.5e, %.5e)", "(zlo, zhi)", ctx.origin.z, ctx.origin.z + ctx.H.m33) << std::endl;
  lout << onika::format_string("- %-15s = (%.5e, %.5e, %.5e)", "(xy, xz, yz)", ctx.H.m12, ctx.H.m13, ctx.H.m23) << std::endl;

  return true;
}

bool LAMMPSDataParser::parse_header(Context& ctx) {

  std::array<bool, static_cast<size_t>(Field::END)> matched_field = {};

  while (!m_file.eof()) {

    current_line() = trim_spaces(current_line());
    if (current_line().empty() || current_line().front() == '#') {
      current_line() = m_file.get_line();
      continue;
    }
    tokenize(current_line(), m_tokens);

    for (size_t i = 0; i < field_readers().size(); ++i) {

      if (matched_field[i])
        continue;

      // auto& [field, parser] = m_field_parsers[i];
      auto& [field, parser] = field_readers()[i];

      if (parser(m_tokens, ctx)) {
        matched_field[i] = true;
        break;
      }
    }

    // if any section keyword is meet we break from parsing headers
    // and go to the next section
    if (matched_field[6]) {
      m_current_section = get_section(m_tokens[0]);
      break;
    }
    current_line() = m_file.get_line();
  }
  return true;
}

bool LAMMPSDataParser::parse_atoms(Context& ctx) {

  // Extract atom style from section header
  tokenize(section_line(), m_tokens);
  size_t index;
  if (match_token_set_any(m_tokens, tok_atom_style, index)) {
    if (m_tokens[index] == "atomic") {
      m_atom_style = AtomStyle::ATOMIC;
    } else {
      m_atom_style = AtomStyle::NONE;
    }
  } else {
    m_atom_style = AtomStyle::NONE;
  }

  // Skip empty/comment lines until first atom data line
  while (current_line().empty() || current_line().front() == '#') {
    if (m_file.eof()) {
      lerr() << "Unexpected EOF while searching for atom data." << std::endl;
      return false;
    }
    current_line() = trim_spaces(m_file.get_line());
  }

  // Tokenize first data line and validate
  tokenize(current_line(), m_tokens);

  if (m_atom_style == AtomStyle::NONE) {
    // TODO: try to guess atom style from
  } else {
    // TODO: check that number of arguments are valid for the given style
  }

  if (m_atom_style == AtomStyle::NONE) {
    lerr() << "Could not determine atom style." << std::endl;
    std::abort();
    return false;
  }

  size_t particle_count = 0;
  AtomLineReader reader = atom_line_readers().at(m_atom_style);

  // pre-allocate particle data
  // ctx.particle_data.reserve(ctx.n_particles);
  ctx.particle_data.resize(ctx.n_particles);

  while (!m_file.eof() && particle_count < ctx.n_particles) {

    current_line() = trim_spaces(current_line());
    
    if (current_line().empty() || current_line().front() == '#') {
      current_line() = m_file.get_line();
      continue;
    }

    tokenize(current_line(), m_tokens);

    if (!reader(m_tokens, ctx, particle_count)) {
      lerr() << "Error parsing atom line at particle #" << particle_count << " : " << current_line()
             << std::endl;
      break;
    }

    ++particle_count;
    current_line() = m_file.get_line();
  }

  if (particle_count != ctx.n_particles) {
    lerr() << "Mismatch in number of particles. Parsed: " << particle_count
           << " Expected: " << ctx.n_particles << std::endl;
    return false;
  }

  m_current_section = Section::IGNORED;
  return true;
}


// TODO:
bool LAMMPSDataParser::parse_masses() {
  m_current_section = Section::IGNORED;
  return true;
}

// TODO:
bool LAMMPSDataParser::parse_velocities(Context&) {
  m_current_section = Section::IGNORED;
  return true;
}

void LAMMPSDataParser::forward_to_next_section() {
  while (!m_file.eof()) {

    current_line() = trim_spaces(current_line());
    if (current_line().empty() || current_line().front() == '#') {
      current_line() = m_file.get_line();
      continue;
    }

    tokenize(current_line(), m_tokens);

    if (match_token_any(m_tokens[0], tok_sections)) {
      m_current_section = get_section(m_tokens[0]);
      section_line() = current_line();
      break;
    }

    current_line() = m_file.get_line();
  }
}

bool read_atom_line_atomic(const TokenSet& tokens, Context& ctx, size_t particle_count) {
  size_t id, type;
  double x, y, z;
  bool ok = tokens_to_num(tokens.at<0>(), id, tokens.at<1>(), type, tokens.at<2>(), x, tokens.at<3>(), y, tokens.at<4>(), z);

  if (ok) {
    ParticleTupleIO& p = ctx.particle_data[particle_count];
    p[field::rx] = x;
    p[field::ry] = y;
    p[field::rz] = z;
    p[field::vx] = 0.0;
    p[field::vy] = 0.0;
    p[field::vz] = 0.0;
    p[field::id] = id;
    p[field::type] = type - 1; // lammps indexing start at one
  }

  return ok;
}

bool read_field_natom(const TokenSet& tokens, Context& ctx) {
  return match_token_any(tokens[1], tok_atoms) && (tokens.size() >= 2) &&
         tokens_to_num(tokens[0], ctx.n_particles);
}

bool read_field_atom_types(const TokenSet& tokens, Context& ctx) {
  if (!match_token_sq(tokens, tok_atom_types))
    return false;
  return tokens_to_num(tokens[0], ctx.n_species);
}

bool read_field_xaxis(const TokenSet& tokens, Context& ctx) {
  double xlo, xhi;
  bool ok = match_token_sq(tokens, tok_xaxis) && tokens_to_num(tokens[0], xlo, tokens[1], xhi);
  if (ok) {
    ctx.H.m11 = xhi - xlo;
    ctx.origin.x = xlo;
  }
  return ok;
}

bool read_field_yaxis(const TokenSet& tokens, Context& ctx) {
  double ylo, yhi;
  bool ok = match_token_sq(tokens, tok_yaxis) && tokens_to_num(tokens[0], ylo, tokens[1], yhi);
  if (ok) {
    ctx.H.m22 = yhi - ylo;
    ctx.origin.y = ylo;
  }
  return ok;
}

bool read_field_zaxis(const TokenSet& tokens, Context& ctx) {
  double zlo, zhi;
  bool ok = match_token_sq(tokens, tok_zaxis) && tokens_to_num(tokens[0], zlo, tokens[1], zhi);
  if (ok) {
    ctx.H.m33 = zhi - zlo;
    ctx.origin.z = zlo;
  }
  return ok;
}

bool read_field_tilt(const TokenSet& tokens, Context& ctx) {
  double xy, xz, yz;
  bool ok = match_token_sq(tokens, tok_tilt) && (tokens.size() >= 6) &&
            tokens_to_num(tokens[0], xy, tokens[1], xz, tokens[2], yz);
  if (ok) {
    ctx.H.m12 = xy;
    ctx.H.m13 = xz;
    ctx.H.m23 = yz;
  }
  return ok;
}

bool read_field_section(const TokenSet& tokens, Context&) { 
  return match_token_set_any(tokens, tok_sections);
}

} // namespace lmpdata

// /* ------------------------------------------------------------------------- */

// namespace lmpdump {
// using namespace re;

// int64_t LAMMPSDumpParser::next() {
//   if (size_t cursor = m_file.tell(); cursor == 0) {
//     m_file.get_line();
//     return static_cast<int64_t>(cursor);
//   } else {
//     return -1;
//   }
// }

// bool LAMMPSDumpParser::parse(Context& ctx) {

//   onika::lout << "Parsing LAMMPS Dump file:" << std::endl;

//   AtomFields atom_fields;

//   // forward to the first valid 'ITEM: NUMBER OF ATOMS'
//   while ((!m_file.eof() &&
//           !std::regex_search(m_current_line.begin(), m_current_line.end(), m_matches[0], re_item_natom))) {
//     m_current_line = m_file.get_line();
//   }

//   // expect to have the number of atom in the next line
//   m_current_line = m_file.get_line();
//   if (!std::regex_search(m_current_line.begin(), m_current_line.end(), m_matches[0], re_uint))
//     return false;
//   ctx.n_particles = static_cast<size_t>(std::stoi(m_matches[0][1].str()));
//   onika::lout << onika::format_string(" - %-15s = %ld", "particle_count", ctx.n_particles) << std::endl;

//   // 'ITEM: BOX BOUNDS' is expected just after NUMBER OF ATOMS
//   m_current_line = m_file.get_line();
//   if (!std::regex_search(m_current_line.begin(), m_current_line.end(), m_matches[0], re_item_bbox))
//     return false;

//   // parsing will be different for orthogonal vs triclinic
//   double xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz;
//   if (std::regex_search(m_current_line.begin(), m_current_line.end(), m_matches[1], re_item_tilt)) {

//     atom_fields.triclinic = true;

//     for (size_t i = 0; i < 3; ++i) {
//       m_current_line = m_file.get_line();
//       if (!std::regex_search(m_current_line.begin(), m_current_line.end(), m_matches[1 + i], re_vec3d))
//         return false;
//     }

//     xlo = std::stod(m_matches[1][1].str());
//     xhi = std::stod(m_matches[1][2].str());

//     ylo = std::stod(m_matches[2][1].str());
//     yhi = std::stod(m_matches[2][2].str());

//     zlo = std::stod(m_matches[3][1].str());
//     zhi = std::stod(m_matches[3][2].str());

//     xy = std::stod(m_matches[1][3].str());
//     xz = std::stod(m_matches[2][3].str());
//     yz = std::stod(m_matches[3][3].str());

//   } else {

//     for (size_t i = 0; i < 3; ++i) {
//       m_current_line = m_file.get_line();
//       if (!std::regex_search(m_current_line.begin(), m_current_line.end(), m_matches[1 + i], re_vec2d))
//         return false;
//     }

//     xlo = std::stod(m_matches[1][1].str());
//     xhi = std::stod(m_matches[1][2].str());

//     ylo = std::stod(m_matches[2][1].str());
//     yhi = std::stod(m_matches[2][2].str());

//     zlo = std::stod(m_matches[3][1].str());
//     zhi = std::stod(m_matches[3][2].str());

//     xy = 0.0;
//     xz = 0.0;
//     yz = 0.0;
//   }

//   onika::lout << onika::format_string(" - %-15s = ", "triclinic") << atom_fields.triclinic << std::endl;

//   xlo -= std::min(xy + xz, std::min(xz, std::min(0., xy)));
//   xhi -= std::max(xy + xz, std::max(xz, std::max(0., xy)));
//   ylo -= std::min(0., yz);
//   yhi -= std::max(0., yz);

//   // vector a
//   ctx.H.m11 = (xhi - xlo);
//   ctx.H.m21 = 0.0;
//   ctx.H.m31 = 0.0;
//   // vector b
//   ctx.H.m12 = xy;
//   ctx.H.m22 = (yhi - ylo);
//   ctx.H.m32 = 0.0;
//   // vector c
//   ctx.H.m13 = xz;
//   ctx.H.m23 = yz;
//   ctx.H.m33 = (zhi - zlo);

//   ctx.origin.x = xlo;
//   ctx.origin.y = ylo;
//   ctx.origin.z = zlo;

//   {
//     // 'ITEM: ATOMS' is expected just after BOX BOUNDS
//     m_current_line = m_file.get_line();
//     if (!std::regex_search(m_current_line.begin(), m_current_line.end(), m_matches[0], re_item_atoms))
//       return false;

//     // parse the properties.
//     // ensure at least x(su) y(su) and z(su) are defined.
//     // if type is not defined we assume type = 0
//     // if id is not defined assume atoms are sorted
//     std::string_view buffer(m_matches[0][0].second, m_current_line.end());

//     // if fields processing fails return
//     if (!process_fields(atom_fields, buffer))
//       return false;

//     onika::lout << onika::format_string(" - %-15s = ", "scaled") << atom_fields.scaled << std::endl;
//     onika::lout << onika::format_string(" - %-15s = ", "unwrapped") << atom_fields.unwrapped << std::endl;

//     // build line regex
//     std::regex re_line_reader = std::regex(build_line_regex_pattern(atom_fields.size));

//     size_t particle_count = 0;
//     ctx.particle_data.reserve(ctx.n_particles);

//     while (!m_file.eof()) {
//       m_current_line = m_file.get_line();

//       if (!std::regex_search(m_current_line.begin(), m_current_line.end(), m_matches[0], re_line_reader))
//         break;

//       size_t type = 0;
//       if (atom_fields.has_type) {
//         size_t index = atom_fields.indices.find(lmpdump::Field::TYPE)->second;
//         type = static_cast<size_t>(std::stod(m_matches[0][1 + index].str())) - 1;
//       }

//       size_t atom_id = particle_count;
//       if (atom_fields.has_id) {
//         size_t index = atom_fields.indices.find(lmpdump::Field::ID)->second;
//         atom_id = static_cast<size_t>(std::stod(m_matches[0][1 + index].str())) - 1;
//       }

//       double x = std::stod(m_matches[0][1 + atom_fields.indices.find(atom_fields.xfield())->second].str());
//       double y = std::stod(m_matches[0][1 + atom_fields.indices.find(atom_fields.yfield())->second].str());
//       double z = std::stod(m_matches[0][1 + atom_fields.indices.find(atom_fields.zfield())->second].str());

//       xyz_fields(x, y, z, atom_fields, ctx);

//       ctx.particle_data.push_back(ParticleTupleIO(x, y, z, atom_id, type));

//       ++particle_count;
//       if (particle_count == ctx.n_particles)
//         break;
//     }

//     if (particle_count != ctx.n_particles)
//       return false;
//   }

//   return true;
// }

// void LAMMPSDumpParser::xyz_fields(double& x, double& y, double& z, const AtomFields& fields, const Context& ctx) {

//   if (!fields.scaled) {
//     return;

//   } else if (!fields.triclinic) {
//     x = x * ctx.H.m11 + ctx.origin.x;
//     y = y * ctx.H.m22 + ctx.origin.y;
//     z = z * ctx.H.m33 + ctx.origin.z;

//   } else {
//     // x = lx * x + xy * y + xz * z + xlo
//     // y = ly * y + xz * z + ylo
//     // z = lz * z + zlo;
//     x = x * ctx.H.m11 + ctx.H.m12 * y + ctx.H.m13 * z + ctx.origin.x;
//     y = y * ctx.H.m22 + ctx.H.m13 * z + ctx.origin.y;
//     z = z * ctx.H.m33 + ctx.origin.z;
//   }
// }

// bool LAMMPSDumpParser::process_fields(AtomFields& fields, const std::string_view& buffer) {
//   sv_regex_iterator it_begin(buffer.begin(), buffer.end(), re_item_fields);
//   sv_regex_iterator it_end = sv_regex_iterator();

//   FieldMap::const_iterator ifield;

//   size_t index = 0;
//   for (; it_begin != it_end; ++it_begin, ++index) {
//     m_matches[1] = *it_begin;

//     ifield = field_map().find(m_matches[1][0].str());

//     if (ifield != field_map().end()) {
//       fields.indices.insert({ifield->second, index});
//     }

//     fields.size++;
//   }

//   // check if type and id field are defined
//   fields.has_type = (fields.indices.find(lmpdump::Field::TYPE) != fields.indices.end());
//   fields.has_type = (fields.indices.find(lmpdump::Field::ID) != fields.indices.end());

//   // check the how the position are given
//   for (size_t i = 0; i < x_priority_groups.size(); ++i) {
//     const auto& group = x_priority_groups[i];
//     if (std::all_of(group.first.begin(), group.first.end(),
//                     [&](const lmpdump::Field& f) { return fields.indices.find(f) != fields.indices.end(); })) {
//       fields.x = group.second;
//       fields.xfields = group.first;

//       if ((fields.x & (X::XS | X::XSU)) != X::NONE)
//         fields.scaled = true;

//       if ((fields.x & (X::XU | X::XSU)) != X::NONE)
//         fields.unwrapped = true;

//       break; // stop at the first valid group
//     }
//   }

//   if (fields.x == X::NONE) {
//     lerr() << "No valid positions field founds in 'ITEM: ATOM': " << buffer << std::endl;
//     return false;
//   }

//   // everything went fine x)
//   return true;
// }

// } // namespace lmpdump

// /* ------------------------------------------------------------------------- */

namespace xyz {

//TODO: Change this
int64_t ExtendedXYZParser::next() {
  if (size_t cursor = m_file.tell(); cursor == 0) {
    m_file.get_line();
    return static_cast<int64_t>(cursor);
  } else {
    return -1;
  }
}

bool ExtendedXYZParser::parse(Context& ctx) {
  onika::lout << "Parsing Extended XYZ file (v2):" << std::endl;

  // xyz files should start with the number of atom but might have a comment at the top
  // so we search for the first non commented line
  while (current_line().empty() || current_line().front() == '#') {
    if (m_file.eof()) {
      lerr() << "Unexpected EOF while searching for atom data." << std::endl;
      return false;
    }
    current_line() = trim_spaces(m_file.get_line());
  }

  // check if the first non comented line is a valid number
  tokenize(current_line(), m_tokens);
  if (!tokens_to_num(m_tokens[0], ctx.n_particles))
    return false;

  // next line should be the comment line with lattice and properties definitions
  current_line() = trim_spaces(m_file.get_line());
  if (!read_exyz_comment_line(current_line(), m_tokens, m_ptokens))
    return false;

  if (!(parse_exyz_lattice(m_tokens, ctx) && parse_exyz_properties(m_ptokens, m_properties)))
    return false;

  // go to atoms lines
  current_line() = m_file.get_line();

  // Skip empty/comment lines until first atom data line
  while (current_line().empty() || current_line().front() == '#') {
    if (m_file.eof()) {
      lerr() << "Unexpected EOF while searching for atom data." << std::endl;
      return false;
    }
    current_line() = trim_spaces(m_file.get_line());
  }

  size_t particle_count = 0;
  ctx.particle_data.resize(ctx.n_particles);

  while (!m_file.eof() && particle_count < ctx.n_particles) {

    current_line() = trim_spaces(current_line());
    
    if (current_line().empty() || current_line().front() == '#') {
      current_line() = m_file.get_line();
      continue;
    }

    tokenize(current_line(), m_tokens);

    if (!parse_line(ctx, m_properties, m_tokens, particle_count)) {
      lerr() << "Error parsing atom line at particle #" << particle_count << " : " << current_line() << std::endl;
      break;
    }

    ++particle_count;
    current_line() = m_file.get_line();
  }

  if (particle_count != ctx.n_particles) {
    lerr() << "Mismatch in number of particles. Parsed: " << particle_count
           << " Expected: " << ctx.n_particles << std::endl;
    return false;
  }

  // print some info after parsing
  lout << onika::format_string("- %-15s = %ld", "particle_count", ctx.n_particles) << std::endl;
  lout << onika::format_string("- %-15s = %ld", "atom_types", ctx.n_species) << std::endl;
  lout << onika::format_string("- %-15s = (%.5e, %.5e)", "(xlo, xhi)", ctx.origin.x, ctx.origin.x + ctx.H.m11) << std::endl;
  lout << onika::format_string("- %-15s = (%.5e, %.5e)", "(ylo, yhi)", ctx.origin.y, ctx.origin.y + ctx.H.m22) << std::endl;
  lout << onika::format_string("- %-15s = (%.5e, %.5e)", "(zlo, zhi)", ctx.origin.z, ctx.origin.z + ctx.H.m33) << std::endl;
  lout << onika::format_string("- %-15s = (%.5e, %.5e, %.5e)", "(xy, xz, yz)", ctx.H.m12, ctx.H.m13, ctx.H.m23) << std::endl;
  return true;
}

} // namespace xyz

/* ------------------------------------------------------------------------- */

std::string remove_leading_chars(const std::string& str, const char c) {
  size_t start = 0;
  while (start < str.size() && str[start] == c) {
    ++start;
  }
  return str.substr(start);
}

const FileInfos guess_file_infos(const std::string& filepath) {

  std::string extension;
  std::filesystem::path path(filepath);

  // remove leading dots return by extension()
  extension = remove_leading_chars(path.extension(), '.');

  // if no extension can't do anything here
  if (extension.empty()) {
    return {.compression = FileCompression::NONE, .format = parser_factory().from_alias("")};
  }

  FileCompression compression = convert_char_to_compression(extension);

  // if first extension is related to compression, get the second one
  if (compression != FileCompression::NONE) {
    extension = remove_leading_chars(path.stem().extension(), '.');
  }

  // return file infos either if format is invalid, this is handled after
  return {.compression = compression, .format = parser_factory().from_alias(extension)};
}

const FileInfos get_file_infos(const std::string& filepath, std::string user_format, std::string user_compression) {

  // never trust the user input
  FileInfos guessed_infos = guess_file_infos(filepath);

  RegisteredFormat format = parser_factory().from_alias(user_format);

  if (!format.is_valid) {
    if (!user_format.empty())
      lout << onika::format_string("[WARN] format='%s' is not a valid format alias\n", user_format.c_str());
    if (!guessed_infos.format.is_valid) {
      lout << onika::format_string("[WARN] Could not deduced format from file extenson\n");
      lout << parser_factory().infos() << std::endl;
      std::abort();
    }
  }

  if (format.is_valid && guessed_infos.format.is_valid)
    lout << onika::format_string("[WARN] Guessed format (%s) does not match provided format (%s)\n",
                                 guessed_infos.format.metadata.name, format.metadata.name);

  FileInfos infos{.compression = user_compression.empty() ? guessed_infos.compression
                                                          : convert_char_to_compression(user_compression),
                  .format = format.is_valid ? format : guessed_infos.format};

  lout << onika::format_string("File infos:\n", filepath.c_str());
  lout << onika::format_string(" - Input file   = '%s' \n", filepath.c_str());
  lout << onika::format_string(" - Atoms format = '%s' \n", infos.format.metadata.name.c_str());
  lout << onika::format_string(" - Compression  = '%s' \n", convert_compression_to_char(infos.compression).c_str());

  return infos;
}

} // namespace ioextra
} // namespace exaStamp
