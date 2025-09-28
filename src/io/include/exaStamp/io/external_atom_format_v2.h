/* -*- c++ -*- ----------------------------------------------------------
   Atoms I/O
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Eliott T. Dubois
------------------------------------------------------------------------- */
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <charconv>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <iostream>
#include <span>
#include <string>
#include <string_view>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

#define AIO_TOKEN_SET_SIZE 128
#define AIO_SPLIT_MAX_COUNT 1000

#define NULL_CHAR '\0'

inline std::string strf(const char* format, ...) {
  char buffer[1024];
  va_list args;
  va_start(args, format);
  vsprintf(buffer, format, args);
  return std::string(buffer);
}

// ------------------------------------------------------------------------- //
#define AIO_POLICY_EXASTAMP_HEADER "onika/math/basic_types.h"
#define AIO_POLICY_HAS_EXASTAMP 0
#ifdef AIO_POLICY_EXASTAMP_HEADER
#if __has_include(AIO_POLICY_EXASTAMP_HEADER)
#undef AIO_POLICY_HAS_EXASTAMP
#define AIO_POLICY_HAS_EXASTAMP 1
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <onika/file_utils.h>
#include <onika/log.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/string_utils.h>

using namespace exanb;

#endif
#endif

#ifdef USE_ZLIB
#include <zconf.h>
#include <zlib.h>
constexpr bool has_zlib = true;
#else
constexpr bool has_zlib = false;
#endif

#ifdef USE_BZIP2
#include <bzlib.h>
constexpr bool has_bzip2 = true;
#else
constexpr bool has_bzip2 = false;
#endif

#ifdef USE_LZMA
#include <lzma.h>
constexpr bool has_lzma = true;
#else
constexpr bool has_lzma = false;
#endif

// ------------------------------------------------------------------------- //

namespace cmp {

// implementation from c++ reference to be use without c++20

template <class T, class U> constexpr bool cmp_equal_int(T t, U u) noexcept {
  if constexpr (std::is_signed_v<T> == std::is_signed_v<U>) {
    return t == u;
  } else if constexpr (std::is_signed_v<T>) {
    // T is signed, U is unsigned
    return t >= 0 && static_cast<std::make_unsigned_t<T>>(t) == u;
  } else {
    // T is signed, U is unsigned
    return u >= 0 && t == static_cast<std::make_unsigned_t<U>>(u);
  }
}

template <class T, class U> constexpr bool cmp_not_equal_int(T t, U u) noexcept { return !cmp_equal_int(t, u); }

template <class T, class U> constexpr bool cmp_less_int(T t, U u) noexcept {
  if constexpr (std::is_signed_v<T> == std::is_signed_v<U>) {
    return t < u;
  } else if constexpr (std::is_signed_v<T>) {
    return t < 0 || static_cast<std::make_unsigned_t<T>>(t) < u;
  } else {
    return u >= 0 && t < static_cast<std::make_unsigned_t<U>>(u);
  }
}

template <class T, class U> constexpr bool cmp_greater_int(T t, U u) noexcept { return cmp_less_int(u, t); }

template <class T, class U> constexpr bool cmp_less_equal_int(T t, U u) noexcept { return !cmp_less_int(u, t); }

template <class T, class U> constexpr bool cmp_greater_equal_int(T t, U u) noexcept { return !cmp_less_int(t, u); }

// wrapper functions

template <class T, class U> constexpr inline bool eq(T t, U u) noexcept { return cmp_equal_int(t, u); }
template <class T, class U> constexpr inline bool ne(T t, U u) noexcept { return cmp_not_equal_int(t, u); }
template <class T, class U> constexpr inline bool lt(T t, U u) noexcept { return cmp_less_int(t, u); }
template <class T, class U> constexpr inline bool le(T t, U u) noexcept { return cmp_less_equal_int(t, u); }
template <class T, class U> constexpr inline bool gt(T t, U u) noexcept { return cmp_greater_int(t, u); }
template <class T, class U> constexpr inline bool ge(T t, U u) noexcept { return cmp_greater_equal_int(t, u); }

template <class T, class U> constexpr auto min(const T& a, const U& b) { return lt(b, a) ? b : a; }

template <class T, class U> constexpr auto max(const T& a, const U& b) { return lt(b, a) ? a : b; }

template <class T, class U>
constexpr inline bool nearly_eq(T t, U u,
                                typename std::common_type<T, U>::type epsilon =
                                    std::numeric_limits<typename std::common_type<T, U>::type>::epsilon()) noexcept {

  using common_t = typename std::common_type<T, U>::type;
  static_assert(std::is_floating_point_v<common_t>, "nearly_equal requires at least one floating-point type");

  common_t a = static_cast<common_t>(t);
  common_t b = static_cast<common_t>(u);
  common_t diff = std::fabs(a - b);
  common_t norm = std::max(std::fabs(a), std::fabs(b));
  return diff <= epsilon * norm;
}

template <class T, class U>
constexpr inline bool nearly_lt(T t, U u,
                                typename std::common_type<T, U>::type epsilon =
                                    std::numeric_limits<typename std::common_type<T, U>::type>::epsilon()) noexcept {

  using common_t = typename std::common_type<T, U>::type;
  return !nearly_eq(t, u, epsilon) && static_cast<common_t>(t) < static_cast<common_t>(u);
}

template <class T, class U>
constexpr inline bool nearly_gt(T t, U u,
                                typename std::common_type<T, U>::type epsilon =
                                    std::numeric_limits<typename std::common_type<T, U>::type>::epsilon()) noexcept {
  return nearly_lt(u, t, epsilon);
}

template <class T, class U>
constexpr inline bool nearly_le(T t, U u,
                                typename std::common_type<T, U>::type epsilon =
                                    std::numeric_limits<typename std::common_type<T, U>::type>::epsilon()) noexcept {
  return !nearly_lt(t, u, epsilon);
}

template <class T, class U>
constexpr inline bool nearly_ge(T t, U u,
                                typename std::common_type<T, U>::type epsilon =
                                    std::numeric_limits<typename std::common_type<T, U>::type>::epsilon()) noexcept {
  return !nearly_lt(u, t, epsilon);
}

template <typename T> constexpr T floor(T x) {
  static_assert(std::is_floating_point_v<T>, "constexpr_floor: T must be floating-point");
  long long xi = static_cast<long long>(x); // truncates toward zero
  return (x < static_cast<T>(xi)) ? static_cast<T>(xi - 1) : static_cast<T>(xi);
}

template <typename T> constexpr T wrap(T x, T max) {
  static_assert(std::is_floating_point_v<T>, "constexpr_wrap: T must be floating-point");
  T r = x - max * cmp::floor(x / max);
  return (r < static_cast<T>(0)) ? r + max : r;
}

} // namespace cmp

namespace enum_traits {

template <typename T>
concept IsEnum = std::is_enum_v<T>;

template <typename T>
concept EnumWithMemberCount = requires { T::Count; } && IsEnum<T>;

template <typename T>
concept EnumWithMemberNone = requires { T::None; } && IsEnum<T>;

template <typename T>
concept EnumWithCountAndNone = EnumWithMemberCount<T> && EnumWithMemberNone<T>;

template <typename T, typename C> inline constexpr T enum_to_type(C EnumMember) { return static_cast<T>(EnumMember); }

template <EnumWithMemberCount T> constexpr size_t enum_size() { return static_cast<size_t>(T::Count); }

template <typename T> constexpr size_t enum_to_index(T t) { return static_cast<size_t>(t); }

template <typename T> constexpr T enum_from_index(size_t index) { return static_cast<T>(index); }

// Helper to convert enum class to underlying type
template <typename Enum> constexpr auto to_underlying(Enum e) noexcept {
  return static_cast<std::underlying_type_t<Enum>>(e);
}

// Generic bitwise operators for any enum class
template <typename Enum> constexpr auto operator|(Enum lhs, Enum rhs) noexcept {
  static_assert(std::is_enum_v<Enum>, "Enum required");
  return static_cast<Enum>(to_underlying(lhs) | to_underlying(rhs));
}

template <typename Enum> constexpr auto operator&(Enum lhs, Enum rhs) noexcept {
  static_assert(std::is_enum_v<Enum>, "Enum required");
  return static_cast<Enum>(to_underlying(lhs) & to_underlying(rhs));
}

template <typename Enum> constexpr auto& operator|=(Enum& lhs, Enum rhs) noexcept {
  lhs = lhs | rhs;
  return lhs;
}

template <typename Enum> constexpr auto& operator&=(Enum& lhs, Enum rhs) noexcept {
  lhs = lhs & rhs;
  return lhs;
}
} // namespace enum_traits

namespace functor {

template <typename Functor>
concept DispatchableFunctor = requires {
  { Functor::N } -> std::convertible_to<size_t>;
};

template <typename Functor>
  requires DispatchableFunctor<Functor>
struct BoolDispatcher {

  static constexpr size_t N = Functor::N;

  template <size_t... Is>
  static constexpr auto deduce_sig_impl(std::index_sequence<Is...>)
      -> decltype(&Functor::template operator()<((void)Is, false)...>);

  struct DeduceSignature {
    using type = decltype(deduce_sig_impl(std::make_index_sequence<N>{}));
  };

  using FnPtr = typename DeduceSignature::type;

  template <size_t I, size_t... Bits> static constexpr FnPtr make_one(std::index_sequence<Bits...>) {
    return &Functor::template operator()<static_cast<bool>((I >> (N - 1 - Bits)) & 1)...>;
  }

  template <size_t... Is> static constexpr auto make(std::index_sequence<Is...>) {
    return std::array<FnPtr, sizeof...(Is)>{make_one<Is>(std::make_index_sequence<N>{})...};
  }

  static constexpr auto table = make(std::make_index_sequence<1u << N>{});

  static size_t index_from_flag(const std::array<bool, N>& flags) {
    size_t idx = 0;
    for (size_t i = 0; i < N; ++i)
      idx = (idx << 1) | flags[i];
    return idx;
  };

  struct Caller {
    Functor functor;
    FnPtr fn;

    template <typename... Args>
    auto operator()(Args&&... args) -> decltype((functor.*fn)(std::forward<Args>(args)...)) {
      return (functor.*fn)(std::forward<Args>(args)...);
    }
  };

  template <typename... Bools> static Caller get(Bools... flags) {
    static_assert(sizeof...(flags) == N);
    return {Functor{}, table[index_from_flag({static_cast<bool>(flags)...})]};
  }
};

} // namespace functor

namespace tokens {

using StringView = std::string_view;
using StringViewSpan = std::span<const StringView>;

/* ------------------------------------------------------------------------- */

inline bool startswith(StringView str, StringView prefix) {
  return (str.size() >= prefix.size()) && (std::memcmp(str.data(), prefix.data(), prefix.size()) == 0);
}

inline bool endswith(StringView str, StringView suffix) {
  const size_t delta = str.size() - suffix.size();
  return (str.size() >= suffix.size()) && (std::memcmp(str.data() + delta, suffix.data(), suffix.size()) == 0);
}

/* ------------------------------------------------------------------------- */

template <size_t Size> struct TokenSetTmpl {

  std::array<StringView, Size> tokens{};
  size_t len = 0;
  StringView* ptr = tokens.data();

  inline void clear() noexcept { len = 0; }
  inline void flush() noexcept {
    std::memset(tokens.data(), 0, len);
    len = 0;
  }

  inline StringView* data() noexcept { return tokens.data(); }
  inline const StringView* data() const noexcept { return tokens.data(); }

  inline StringView operator[](size_t i) noexcept { return tokens[i]; }
  inline const StringView operator[](size_t i) const noexcept { return tokens[i]; }

  template <size_t index> inline constexpr StringView& at() noexcept {
    static_assert(index < Size);
    return ptr[index];
  }

  template <size_t index> inline constexpr StringView& at() const noexcept {
    static_assert(index < Size);
    return ptr[index];
  }

  inline constexpr size_t size() const noexcept { return len; }
  inline constexpr size_t capacity() const noexcept { return Size; }
  inline constexpr bool full() const noexcept { return size() >= capacity(); }
  inline constexpr bool atleast(size_t n) const noexcept { return size() >= n; }

  inline void push_back(StringView token) noexcept {
    if (!full()) {
      tokens[len++] = token;
    }
  }

  inline void push_back_unchecked(StringView token) noexcept { tokens[len++] = token; }

  StringView* begin() noexcept { return tokens.begin(); }
  const StringView* begin() const noexcept { return tokens.begin(); }
  StringView* end() noexcept { return tokens.begin() + len; }
  const StringView* end() const noexcept { return tokens.begin() + len; }
};

using TokenSet = TokenSetTmpl<AIO_TOKEN_SET_SIZE>;

template <size_t N> using TokenNeedles = std::array<StringView, N>;

/* ------------------------------------------------------------------------- */

template <char... Chars> struct IsChar {
  constexpr bool operator()(char c) const noexcept { return ((c == Chars) || ...); }
  constexpr bool operator()(const char* cptr) const noexcept { return cptr && ((*cptr == Chars) || ...); }
  // constexpr bool operator()(const char* cptr) const noexcept { return cptr && this(*cptr); } //
};

struct IsSpace : IsChar<' ', '\t', '\n', '\r'> {};
struct SpaceDelimiter : IsChar<' ', '\t'> {};
struct CommaDelimiter : IsChar<','> {};
struct ColumnDelimiter : IsChar<':'> {};
struct SemiColumnDelimiter : IsChar<';'> {};
struct DoubleQuoteDelimiter : IsChar<'"'> {};
struct SingleQuoteDelimiter : IsChar<'\''> {};
struct EqualSignDelimiter : IsChar<'='> {};
struct IsHashComment : IsChar<'#'> {};

template <typename Predicate> inline StringView ltrim(StringView str, Predicate&& predicate) noexcept {
  while (!str.empty() && predicate(str.front())) {
    str.remove_prefix(1);
  }
  return str;
}

template <typename Predicate> inline StringView rtrim(StringView str, Predicate&& predicate) noexcept {
  while (!str.empty() && predicate(str.back())) {
    str.remove_suffix(1);
  }
  return str;
}

template <typename Predicate> inline StringView trim(StringView str, Predicate&& predicate) noexcept {
  return rtrim(ltrim(str, std::forward<Predicate>(predicate)), std::forward<Predicate>(predicate));
}

inline StringView trim_spaces(StringView str) noexcept { return trim(str, IsSpace{}); }

template <typename Predicate> inline StringView cut_at_first(StringView str, Predicate&& predicate) noexcept {
  const char* begin = str.data();
  const char* end = begin + str.size();
  const char* it = begin;
  while (it < end && !predicate(it))
    ++it;
  return StringView{begin, static_cast<size_t>(it - begin)};
}

template <typename Predicate> inline bool skipline_tmpl(StringView str, Predicate&& predicate) noexcept {
  return str.empty() || predicate(str.front());
}

inline bool skipline(StringView str) noexcept { return skipline_tmpl(str, IsHashComment{}); }

template <typename Delimiter> struct needs_trim : std::true_type {};

template <char... Chars> struct needs_trim<IsChar<Chars...>> {
  static constexpr bool value = !((IsSpace{}(Chars)) || ...);
};

template <typename Callback, typename Delimiter>
inline void split(StringView str, Delimiter&& delimiter, Callback&& callback,
                  size_t max_count = AIO_SPLIT_MAX_COUNT) noexcept {
  size_t count = 0;
  const char* data = str.data();
  const char* end = data + str.size();
  while (data < end && cmp::lt(count, max_count)) {
    while (data < end && delimiter(data))
      ++data;
    const char* token_start = data;
    while (data < end && !delimiter(data))
      ++data;
    if (token_start != data) {
      StringView tok{token_start, static_cast<size_t>(data - token_start)};
      if constexpr (needs_trim<std::decay_t<Delimiter>>::value)
        tok = trim_spaces(tok);
      callback(tok);
      ++count;
    }
  }
}

template <typename Delimiter, size_t N>
inline void tokenize_impl(const StringView str, TokenSetTmpl<N>& tokens, Delimiter&& delimiter,
                          size_t max_count = N) noexcept {
  tokens.clear();
  split(str, delimiter, [&tokens](StringView token) { tokens.push_back_unchecked(token); }, cmp::min(max_count, N));
}

inline void tokenize(const StringView str, TokenSet& tokens) { return tokenize_impl(str, tokens, SpaceDelimiter{}); }

inline void tokenize(const StringView str, TokenSet& tokens, size_t max_count) {
  return tokenize_impl(str, tokens, SpaceDelimiter{}, max_count);
}

/* ------------------------------------------------------------------------- */

template <bool CaseSensitive> struct CharComparator {
  constexpr inline bool operator()(char a, char b) const noexcept {
    if constexpr (CaseSensitive) {
      return a == b;
    } else {
      if (a >= 'A' && a <= 'Z')
        a += 32;
      if (b >= 'A' && b <= 'Z')
        b += 32;
      return (a == b);
    }
  }
};

template <bool CaseSensitive> struct StringViewComparator {
  static constexpr CharComparator<false> comparator{};
  constexpr inline bool operator()(const StringView a, const StringView b) {
    if (a.size() != b.size())
      return false;
    if constexpr (CaseSensitive) {
      return std::memcmp(a.data(), b.data(), a.size()) == 0;
    } else {
      for (size_t i = 0; i < a.size(); ++i) {
        if (!comparator(a[i], b[i]))
          return false;
      }
      return true;
    }
  }
};

template <bool CaseSensitive>
inline bool token_match_impl(const StringView token, const StringView needle,
                             StringViewComparator<CaseSensitive>&& comparator) {
  return comparator(token, needle);
}

inline bool token_match(const StringView token, const StringView needle) {
  return token_match_impl(token, needle, StringViewComparator<true>{});
}

inline bool token_match_no_case(const StringView token, const StringView needle) {
  return token_match_impl(token, needle, StringViewComparator<false>{});
}

// Return true if 'token' match any of the needles
// Stop at the first match and set 'index' to the matching index in 'needles'.
template <bool CaseSensitive>
inline bool token_match_any_impl(const StringView token, const StringViewSpan needles, size_t& index,
                                 StringViewComparator<CaseSensitive>&& comparator) {
  const size_t n = needles.size();

  if (n == 0)
    return false;

  for (index = 0; index < n; ++index) {
    if (comparator(token, needles[index]))
      return true;
  }

  return false;
}

inline bool token_match_any(const StringView token, const StringViewSpan needles, size_t& index) {
  return token_match_any_impl(token, needles, index, StringViewComparator<true>{});
}

inline bool token_match_any(const StringView token, const StringViewSpan needles) {
  size_t index{};
  return token_match_any_impl(token, needles, index, StringViewComparator<true>{});
}

inline bool token_match_any_no_case(const StringView token, const StringViewSpan needles, size_t& index) {
  return token_match_any_impl(token, needles, index, StringViewComparator<false>{});
}

inline bool token_match_any_no_case(const StringView token, const StringViewSpan needles) {
  size_t index{};
  return token_match_any_impl(token, needles, index, StringViewComparator<false>{});
}

// Return true when any token in 'tokens' match any needle in 'needles'
// Stop at the first match and provide i and j indices of the token and needles respectively.
inline bool token_set_match_any(const StringViewSpan tokens, const StringViewSpan needles, size_t& ii, size_t& jj) {
  const size_t token_count = tokens.size();
  const size_t needle_count = needles.size();

  if (needle_count == 0) {
    ii = 0;
    jj = 0;
    return false;
  }

  for (ii = 0; ii < token_count; ++ii) {
    for (jj = 0; jj < needle_count; ++jj) {
      if (tokens[ii] == needles[jj]) {
        return true;
      }
    }
  }

  return false;
}

inline bool token_set_match_any(StringViewSpan tokens, StringViewSpan needles) {
  size_t i, j;
  return token_set_match_any(tokens, needles, i, j);
}

// Returns true if the sequence of needles is present in the sequence of tokens
inline bool token_set_match_seq(const StringViewSpan tokens, const StringViewSpan needles, size_t& index) {
  const size_t token_count = tokens.size();
  const size_t needle_count = needles.size();

  if (needle_count == 0 || token_count < needle_count) {
    return false;
  }

  for (size_t i = 0; i <= token_count - needle_count; ++i) {
    bool match = true;
    for (size_t j = 0; j < needle_count; ++j) {
      if (tokens[i + j] != needles[j]) {
        match = false;
        break;
      }
    }
    if (match) {
      index = i;
      return true;
    }
  }
  return false;
}

inline bool token_set_match_seq(StringViewSpan tokens, StringViewSpan needles) {
  size_t index;
  return token_set_match_seq(tokens, needles, index);
}

template <size_t N> inline bool token_match_any(const StringView token, const TokenNeedles<N>& needles) {
  return token_match_any(token, StringViewSpan(needles.data(), N));
}

template <size_t N> inline bool token_match_any(const StringView token, const TokenNeedles<N>& needles, size_t& index) {
  return token_match_any(token, StringViewSpan(needles.data(), N), index);
}

template <> inline bool token_match_any(const StringView token, const TokenNeedles<1>& needles, size_t& index) {
  index = 0;
  return token_match(token, std::get<0>(needles));
}

template <size_t N, size_t M>
inline bool token_set_match_any(const TokenSetTmpl<M>& tokens, const TokenNeedles<N>& needles) {
  return token_set_match_any(StringViewSpan(tokens.data(), tokens.size()), StringViewSpan(needles.data(), N));
}

template <size_t N, size_t M>
inline bool token_set_match_any(const TokenSetTmpl<M>& tokens, const TokenNeedles<N>& needles, size_t& i, size_t& j) {
  return token_set_match_any(StringViewSpan(tokens.data(), tokens.size()), StringViewSpan(needles.data(), N), i, j);
}

template <size_t N, size_t M>
inline bool token_set_match_seq(const TokenSetTmpl<M>& tokens, const TokenNeedles<N>& needles) {
  return token_set_match_seq(StringViewSpan(tokens.data(), tokens.size()), StringViewSpan(needles.data(), N));
}

template <size_t N, size_t M>
inline bool token_set_match_seq(const TokenSetTmpl<M>& tokens, const TokenNeedles<N>& needles, size_t& index) {
  return token_set_match_seq(StringViewSpan(tokens.data(), tokens.size()), StringViewSpan(needles.data(), N), index);
}

namespace numeric {

template <typename T, typename = void> struct TokenParser;

template <typename T> struct TokenParser<T, std::enable_if_t<std::is_floating_point_v<T>>> {
  static constexpr bool parse(const char* begin, size_t size, T& rval) {
    const char* end = begin + size;
    auto [ptr, err] = std::from_chars(begin, end, static_cast<T&>(rval), std::chars_format::general);
    return (err == std::errc() && ptr == end);
  }
};

template <typename T> struct TokenParser<T, std::enable_if_t<std::is_integral_v<T>>> {

  static constexpr T min = std::numeric_limits<T>::min();
  static constexpr T max = std::numeric_limits<T>::max();

  static constexpr bool parse(const char* begin, size_t size, T& rval) {
    const char* end = begin + size;
    // try to parse as integer directly
    auto [ptr, err] = std::from_chars(begin, end, static_cast<T&>(rval));
    if (err == std::errc() && ptr == end)
      return true;

    // try to parse as double (needed for scientific notation)
    double tmp = 0;
    auto [fptr, fec] = std::from_chars(begin, end, tmp, std::chars_format::general);
    if (fec != std::errc{} || fptr != end)
      return false;
    if (tmp != static_cast<double>(static_cast<T>(tmp)))
      return false;
    if (tmp < static_cast<double>(min) || tmp > static_cast<double>(max))
      return false;
    rval = static_cast<T>(tmp);
    return true;
  }
};

template <typename T> inline bool parse_single_token(const StringView token, T& rval) {
  static_assert(std::is_arithmetic_v<T>, "parse_single_token_impl requires arithmetic type");
  return TokenParser<T>::parse(token.begin(), token.size(), rval);
}

inline bool parse_tokens() noexcept { return true; }

template <typename TokenType, typename ValueType, typename... Args>
inline bool parse_tokens(const TokenType& token, ValueType& value, Args&&... args) noexcept {
  if constexpr (sizeof...(args) == 0) {
    return parse_single_token(token, value);
  } else {
    if (parse_single_token(token, value)) {
      return parse_tokens(std::forward<Args>(args)...);
    }
    return false;
  }
}
} // namespace numeric
} // namespace tokens

namespace properties_traits {
using namespace enum_traits;

template <typename Type, typename Enum, Enum Type::* Member>
concept TypeHasEnumMember = requires(Type t) {
  { t.*Member } -> std::same_as<Enum&>;
};

template <typename Type, IsEnum Enum, Enum Type::* Member>
  requires TypeHasEnumMember<Type, Enum, Member> && EnumWithCountAndNone<Enum>
struct GenericContainer {

  using enum_t = Enum;
  static constexpr size_t capacity = enum_size<Enum>();
  size_t len = 0;
  std::array<Type, capacity> data{};

  template <Enum e> static constexpr void assert_size() {
    static_assert(enum_to_index(e) < capacity, "PropertiesContainer::get : out-of-bound");
  }

  template <Enum e> constexpr const Type& get() const {
    assert_size<e>();
    return data[enum_to_index(e)];
  }

  const Type& get(Enum e) const {
    if (enum_to_index(e) < capacity) {
      return data[enum_to_index(e)];
    } else {
      return data[enum_to_index(Enum::None)];
    }
  }

  template <Enum... args> constexpr bool has() const { return ((get<args>().*Member == args) && ...); }

  template <typename... Args> bool has(Args... args) {
    static_assert((std::is_same_v<Args, Enum> && ...), "Only E enum type allowed");
    return ((get(args).*Member == args) && ...);
  }

  template <Enum... args> constexpr bool any() const { return ((get<args>().*Member == args) || ...); }

  template <typename... Args> constexpr void expect(Args... args) {
    static_assert((std::is_same_v<Args, Enum> && ...), "Only E enum type allowed");
    ((len += (has(args)) ? static_cast<size_t>(1) : static_cast<size_t>(0)), ...);
  }

  inline void set(Type t) {
    if (enum_to_index<Enum>(t.*Member) < capacity) {
      data[enum_to_index<Enum>(t.*Member)] = t;
    }
  }

  inline void clear() {
    len = 0;
    std::fill(data.begin(), data.end(), Type{});
  }

  inline size_t size(void) const { return len; }
  inline size_t size(void) { return len; }
};
} // namespace properties_traits

namespace aio {

namespace details {

inline uint64_t file_size(const char* path) {
  struct stat st;
  if (stat(path, &st) != 0)
    std::abort();
  return static_cast<std::uint64_t>(st.st_size);
}

inline bool file_exists(const char* path) { return std::filesystem::exists(path); }

template <typename T>
  requires std::is_integral_v<T> && std::is_unsigned_v<T>
inline T pack_string(std::string_view s) {
  constexpr size_t MAX_CHARS = sizeof(T);
  T value = 0;
  const char* ptr = s.data();
  for (size_t i = 0; i < MAX_CHARS; ++i) {
    uint8_t c = (i < s.size()) * static_cast<uint8_t>(ptr[i]);
    value = (value << 8) | c; // shift left, append char
  }
  return value;
}

} // namespace details

namespace shared_types {
struct Vec3d {
  double x;
  double y;
  double z;
};

struct Mat3d {
  double m11, m12, m13;
  double m21, m22, m23;
  double m31, m32, m33;
};

struct IJK {
  ssize_t i;
  ssize_t j;
  ssize_t k;
};
} // namespace shared_types

namespace base_types {

template <typename T, size_t DIM> struct IArray {
  virtual ~IArray() = default;
  virtual T& operator()(size_t i, size_t j) = 0;
  virtual const T& operator()(size_t i, size_t j) const = 0;
  virtual T* data() = 0;
  virtual const T* data() const = 0;
  virtual size_t size() const = 0;
  virtual void resize(size_t n) = 0;
};

template <typename T, size_t DIM> class Array : public IArray<T, DIM> {
public:
  inline T& operator()(size_t i, size_t j) override { return m_data.data()[i * DIM + j]; }
  inline const T& operator()(size_t i, size_t j) const override { return m_data.data()[i * DIM + j]; }
  inline T* data() override { return m_data.data(); }
  inline const T* data() const override { return m_data.data(); }
  inline size_t size() const override { return m_data.size(); }
  inline void resize(size_t n) override { m_data.resize(DIM * n); }

private:
  std::vector<T> m_data;
};

template <typename T, size_t DIM> class ArrayView : public IArray<T, DIM> {
public:
  ArrayView(T* ptr, size_t n) : m_ptr(ptr), m_size(n) {}

  inline T& operator()(size_t i, size_t j) override { return m_ptr[i * DIM + j]; }
  inline const T& operator()(size_t i, size_t j) const override { return m_ptr[i * DIM + j]; }
  inline T* data() override { return m_ptr; }
  inline const T* data() const override { return m_ptr; }
  inline size_t size() const override { return m_size; }
  inline void resize(size_t n) override {
    if (n != m_size)
      std::abort();
  }

private:
  T* m_ptr;
  size_t m_size;
};

template <typename T, size_t DIM> class ParticleArray {
public:
  static constexpr size_t D = DIM;

  ParticleArray() : m_self(std::make_unique<Array<T, DIM>>()) {}

  template <typename ArrayType>
    requires std::derived_from<ArrayType, IArray<T, DIM>>
  ParticleArray(ArrayType array) : m_self(std::make_unique<ArrayType>(std::move(array))) {}

  // move only
  ParticleArray(ParticleArray&&) noexcept = default;
  ParticleArray& operator=(ParticleArray&&) noexcept = default;
  ParticleArray(const ParticleArray&) = delete;
  ParticleArray& operator=(const ParticleArray&) = delete;

  // forward API
  inline T& operator()(size_t i, size_t j) { return (*m_self)(i, j); }
  inline const T& operator()(size_t i, size_t j) const { return (*m_self)(i, j); }
  inline T* data() { return m_self->data(); }
  inline const T* data() const { return m_self->data(); }
  inline size_t size() const { return m_self->size(); }
  inline void resize(size_t n) { m_self->resize(n); }

private:
  std::unique_ptr<IArray<T, DIM>> m_self;
};

template <size_t N> struct SpeciesMapTmpl {
  static constexpr size_t MAX = N;
  size_t count = 0;                      // number of registered type
  ssize_t max_type = -1;                 // maximum index of registered type.
  std::array<uint32_t, N> keys{};        // packed symbol keys
  std::array<ssize_t, N> types{};        // type ID for each index
  std::array<std::string, N> symbols{};  // string symbol
  std::array<size_t, N> type_to_index{}; // map type -> array index

  uint32_t cached_key = 0xffffffff;
  ssize_t cached_type = -1;
  size_t cached_index = -1;

  std::array<size_t, N> counts{};
  std::array<size_t, N> offsets{};
  size_t total_atoms = 0;
  bool offsets_valid = true;

  inline ssize_t next_type() const noexcept { return max_type + 1; }

  // -- Registration --
  void add_species(std::string_view symbol, ssize_t type = -1) {

    if (cmp::gt(count, N) || cmp::gt(type, N))
      return;

    uint32_t key = details::pack_string<uint32_t>(symbol);

    // check for duplicate keys
    for (size_t i = 0; i < count; ++i) {
      if (keys[i] == key)
        return;
    }

    if (cmp::lt(type, 0))
      type = next_type();

    size_t i = count++;
    keys[i] = key;
    types[i] = type;
    symbols[i] = std::string(symbol);
    type_to_index[type] = i;
    max_type = cmp::max(max_type, type);
  }

  // -- Lookup --

  inline ssize_t get_type(std::string_view symbol) {

    ssize_t type = -1;
    uint32_t key = details::pack_string<uint32_t>(symbol);

    if (cmp::eq(key, cached_key))
      return cached_type;

    for (size_t i = 0; i < count; ++i) {
      if (keys[i] == key)
        type = types[i];
    }

    cached_key = key;
    cached_type = type;

    return type;
  }

  inline std::string_view get_symbol(ssize_t type) {
    static constexpr const char* fallback = "Xx";
    if (cmp::lt(type, 0) || cmp::ge(type, N))
      return fallback;
    size_t idx = type_to_index[type];
    if (cmp::ge(idx, count))
      return fallback;
    return symbols[idx];
  }
  // --- Counts and offsets ---
  void set_type_count(ssize_t type, size_t n) {
    if (cmp::lt(type, 0) || cmp::ge(type, N))
      return;
    counts[type] = n;
    offsets_valid = false;
  }

  size_t get_count(ssize_t type) {
    if (cmp::lt(type, 0) || cmp::ge(type, N))
      return 0;
    return counts[type];
  }

  void compute_count_offsets() {
    size_t sum = 0;
    for (ssize_t t = 0; t <= max_type; ++t) {
      offsets[t] = sum;
      sum += counts[t];
    }
    total_atoms = sum;
    offsets_valid = true;
  }

  ssize_t type_from_index(size_t idx) {
    if (!offsets_valid || cmp::ge(idx, total_atoms))
      return -1;

    for (ssize_t t = 0; t <= max_type; ++t) {
      if (cmp::lt(idx, offsets[t] + counts[t]))
        return t;
    }

    return -1; // should never happend
  }

  size_t size() const noexcept { return count; }
  size_t atoms_total() const noexcept { return total_atoms; }
};

template <size_t N> struct SpeciesMapTmplOld {
private:
  struct SumCounter {
    std::array<size_t, N> m_count{};
    std::array<size_t, N> m_offsets{};

    inline void set(size_t type, size_t count) {
      m_count[type] = count;
      compute_offset();
    }

    inline size_t get(size_t type) { return m_count[type]; }

    inline size_t from_index(size_t index) {
      for (size_t type = 1; type < N; ++type) {
        if (index < m_offsets[type])
          return type - 1;
      }
      return 0;
    }

    inline void compute_offset() {
      size_t sum = 0;
      for (size_t type = 0; type < N; ++type) {
        m_offsets[type] = sum;
        sum += m_count[type];
      }
    }
  };

  size_t m_count = 0;

  std::array<uint32_t, N> m_keys{};
  std::array<ssize_t, N> m_type{};
  std::array<size_t, N> m_index{};

  std::array<std::string, N> m_symbols{};

  uint32_t m_cached_key = 0xffffffff;
  int m_cached_type = -1;
  int m_type_mx = -1;

public:
  static constexpr size_t MAX = N;
  SumCounter counter{};

  inline void set(const std::string& symbol, ssize_t type) {
    if (m_count > N || static_cast<size_t>(type) > N)
      throw std::runtime_error("MAX number of species.");

    uint32_t key = details::pack_string<uint32_t>(symbol);
    for (size_t i = 0; i < m_count; ++i)
      if (m_keys[i] == key)
        throw std::runtime_error("Can't register 2 type the same species.");

    size_t index = m_index[type];
    if (index == 0) {
      m_count++;
      index = m_count;
    }

    m_keys[index] = details::pack_string<uint32_t>(symbol);
    m_type[index] = type;
    m_symbols[index] = symbol;
    m_index[type] = index;
    m_type_mx = std::max(m_type_mx, static_cast<int>(type));
  }

  inline void add(const std::string& symbol) { this->set(symbol, m_type_mx + 1); }

  // return the type given the symbol
  inline ssize_t get(std::string_view symbol) {
    uint32_t key = details::pack_string<uint32_t>(symbol);
    if (key == m_cached_key)
      return m_cached_type;
    ssize_t type = -1;
    for (size_t i = 0; i < m_count + 1; ++i) {
      if (m_keys[i] == key) {
        type = m_type[i];
        break;
      }
    }
    m_cached_key = key;
    m_cached_type = type;
    return type;
  };

  // return a symbol given the type
  inline std::string_view get(ssize_t type) {
    size_t idx = m_index[type];
    if (m_keys[idx] == 0)
      return std::string_view("Xx");
    return std::string_view(m_symbols[idx]);
  }

  inline size_t count() { return m_count; }
};

using SpeciesMap = SpeciesMapTmpl<16>;

} // namespace base_types

namespace policy {

enum class Policy {
  None,
  Internal,
  ExaStamp,
};

template <Policy> struct TraitsSelector;

// Traits definition
struct InternalTraits {
  using Vec3d = typename shared_types::Vec3d;
  using Mat3d = typename shared_types::Mat3d;
  using IJK = typename shared_types::IJK;

  struct ContextDataBuffer {
    size_t nat{}, ntypes{};
    Mat3d cell{};
    Vec3d origin{};
    base_types::ParticleArray<double, 3> positions;
    base_types::ParticleArray<double, 3> velocities{};
    base_types::ParticleArray<double, 3> forces{};
    base_types::ParticleArray<int, 1> types{};
  };

  static inline void allocate(ContextDataBuffer& buf, size_t n) {
    buf.positions.resize(n);
    buf.types.resize(n);
    buf.velocities.resize(n);
  }

  static inline void assign_position(ContextDataBuffer& buf, size_t i, double x, double y, double z) {
    double* __restrict p = buf.positions.data();
    p[i * 3 + 0] = x;
    p[i * 3 + 1] = y;
    p[i * 3 + 2] = z;
  }

  static inline void assign_velocity(ContextDataBuffer& buf, size_t i, double vx, double vy, double vz) {
    double* __restrict v = buf.velocities.data();
    v[i * 3 + 0] = vx;
    v[i * 3 + 1] = vy;
    v[i * 3 + 2] = vz;
  }

  static inline void assign_type(ContextDataBuffer& buf, size_t i, size_t type) {
    buf.types.data()[i] = static_cast<int>(type);
  }
};

template <> struct TraitsSelector<Policy::Internal> {
  using type = InternalTraits;
};

#ifdef AIO_POLICY_HAS_EXASTAMP

struct ExaStampTraits {
  using Vec3d = onika::math::Vec3d;
  using Mat3d = onika::math::Mat3d;
  using IJK = onika::math::IJK;
  using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_vx, field::_vy,
                                                   field::_vz, field::_id, field::_type>;
  using ParticleData = std::vector<ParticleTupleIO>;

  struct ContextDataBuffer {
    size_t nat{}, ntypes{};
    Mat3d cell{};
    Vec3d origin{};

    // base_types::ParticleArray<double, 3> positions;
    // base_types::ParticleArray<double, 3> velocities{};
    // base_types::ParticleArray<double, 3> forces{};
    // base_types::ParticleArray<int, 1> types{};
    ParticleData particles;
  };

  static inline void allocate(ContextDataBuffer& buf, size_t n) { buf.particles.resize(n); }

  static inline void assign_position(ContextDataBuffer& buf, size_t i, double x, double y, double z) {
    ParticleTupleIO& p = buf.particles[i];
    p[field::rx] = x;
    p[field::ry] = y;
    p[field::rz] = z;
    p[field::id] = i;
  }

  static inline void assign_velocity(ContextDataBuffer& buf, size_t i, double vx, double vy, double vz) {
    ParticleTupleIO& p = buf.particles[i];
    p[field::vx] = vx;
    p[field::vy] = vy;
    p[field::vz] = vz;
  }

  static inline void assign_type(ContextDataBuffer& buf, size_t i, size_t type) {
    buf.particles[i][field::type] = static_cast<ssize_t>(type);
  }

  static inline void panic(const std::string& msg, const char* file, const char* func, int line) {
    static constexpr const char* str = "ERR: [%s:%s:%d] %s\n";
    const char* fpath = std::filesystem::path(file).filename().c_str();
    onika::lout << strf(str, fpath, func, line, msg.c_str()) << std::flush;
  };

  static inline void log_infos(const std::string& msg) { onika::lout << msg; };
  static inline void log_warn(const std::string& msg) { onika::lout << "[WARN] " << msg; };
};

template <> struct TraitsSelector<Policy::ExaStamp> {
  using type = ExaStampTraits;
};

#endif

// Compile-time choice of the policy
#define AIO_PP_STR_HELPER(x) #x
#define AIO_PP_STR(x) AIO_PP_STR_HELPER(x)

#ifndef AIO_POLICY
#define AIO_POLICY Internal
#endif

constexpr Policy policy_from_string(std::string_view name) {
  if (name == "Internal" || name == "Python")
    return Policy::Internal;
  else if (name == "ExaStamp")
    return Policy::ExaStamp;
  else
    return Policy::None;
}

inline constexpr bool has_exastamp = (AIO_POLICY_HAS_EXASTAMP != 0);
inline constexpr Policy active_policy = policy_from_string(AIO_PP_STR(AIO_POLICY));

static_assert(active_policy != Policy::None, "Invalid AIO_POLICY. Must be one of: Internal, ExaStamp, Lammps, Python.");
static_assert(!(active_policy == Policy::ExaStamp && !has_exastamp),
              "AIO_POLICY=ExaStamp was selected but AIO_POLICY_EXASTAMP_HEADER was not provided or "
              "header is missing.");

using Traits = typename TraitsSelector<active_policy>::type;

#undef AIO_PP_STR_HELPER
#undef AIO_PP_STR

template <typename Policy>
concept HasContextDataBuffer = requires { typename Policy::ContextDataBuffer; };

template <typename Buffer, typename Policy>
concept ContextDataBufferLike = requires(Buffer buf) {
  { buf.nat } -> std::convertible_to<size_t>;
  { buf.ntypes } -> std::convertible_to<size_t>;
  { buf.cell } -> std::convertible_to<typename Policy::Mat3d>;
  { buf.origin } -> std::convertible_to<typename Policy::Vec3d>;
};

template <typename Policy>
concept HasPanic = requires(const std::string& msg, const char* file, const char* func, int line) {
  { Policy::panic(msg, file, func, line) } -> std::same_as<void>;
};

template <typename Policy>
concept HasLogInfo = requires(const std::string& msg) {
  { Policy::log_infos(msg) } -> std::same_as<void>;
};

template <typename Policy>
concept HasLogWarn = requires(const std::string& msg) {
  { Policy::log_warn(msg) } -> std::same_as<void>;
};

template <typename Policy>
concept ValidPolicy = requires {
  typename Policy::Vec3d;
  typename Policy::Mat3d;
  typename Policy::IJK;
  typename Policy::ContextDataBuffer;
} && ContextDataBufferLike<typename Policy::ContextDataBuffer, Policy>;

template <typename Policy>
[[noreturn]] inline void panic(const std::string& msg, const char* file, const char* func, int line) {
  if constexpr (HasPanic<Policy>) {
    Policy::panic(msg, file, func, line);
  } else {
    static constexpr const char* str = "ERR: [%s:%s:%d] %s\n";
    const char* fpath = std::filesystem::path(file).filename().c_str();
    std::cout << strf(str, fpath, func, line, msg.c_str()) << std::flush;
  }
  std::abort();
};

template <typename Policy> inline void log_info(const std::string& msg) {
  if constexpr (HasLogInfo<Policy>) {
    Policy::log_infos(msg);
  } else {
    std::cout << msg;
  }
};

template <typename Policy> inline void log_warn(const std::string& msg) {
  if constexpr (HasLogWarn<Policy>) {
    Policy::log_warn(msg);
  } else {
    std::cout << "[WARN] " << msg;
  }
};

template <Policy P> constexpr void assert_policy_impl() {
  static_assert(active_policy == P, "Active policy do not match requested policy.");
}

} // namespace policy

#define ASSERT_AIO_POLICY(P)                                                                                           \
  static_assert(aio::policy::active_policy == aio::policy::Policy::P, "Active policy do not match requested policy.");

#define PANIC(...)                                                                                                     \
  do {                                                                                                                 \
    aio::policy::panic<aio::policy::Traits>(strf(__VA_ARGS__), __FILE__, __func__, __LINE__);                          \
  } while (0);

#define linfo(...)                                                                                                     \
  do {                                                                                                                 \
    aio::policy::log_info<aio::policy::Traits>(strf(__VA_ARGS__));                                                     \
  } while (0);

#define lwarn(...)                                                                                                     \
  do {                                                                                                                 \
    aio::policy::log_warn<aio::policy::Traits>(strf(__VA_ARGS__));                                                     \
  } while (0);

// add policy-dependent types
namespace base_types {

#define aio_policy_deftype(type) using type = typename policy::Traits::type;

aio_policy_deftype(Vec3d);
aio_policy_deftype(Mat3d);
aio_policy_deftype(IJK);

} // namespace base_types

namespace matrix_traits {

using namespace base_types;

template <typename T>
concept Vec3dLike = requires(T v) {
  sizeof(T) == sizeof(double) * 3;
  std::is_trivially_copyable_v<T>;
  { v.x } -> std::convertible_to<double>;
  { v.y } -> std::convertible_to<double>;
  { v.z } -> std::convertible_to<double>;
};

template <typename T>
concept IJKLike = requires(T i) {
  std::is_trivially_copyable_v<T>;
  sizeof(T) == sizeof(ssize_t) * 3;
  { i.i } -> std::convertible_to<ssize_t>;
  { i.j } -> std::convertible_to<ssize_t>;
  { i.k } -> std::convertible_to<ssize_t>;
};

template <typename T>
concept QuaternionLike = requires(T q) {
  sizeof(T) == sizeof(double) * 4;
  std::is_trivially_copyable_v<T>;
  { q.w } -> std::convertible_to<double>;
  { q.x } -> std::convertible_to<double>;
  { q.y } -> std::convertible_to<double>;
  { q.z } -> std::convertible_to<double>;
};

template <typename T>
concept Mat3dLike = requires(T m) {
  sizeof(T) == sizeof(double) * 9;
  std::is_trivially_copyable_v<T>;
  { m.m11 } -> std::convertible_to<double>;
  { m.m12 } -> std::convertible_to<double>;
  { m.m13 } -> std::convertible_to<double>;
  { m.m21 } -> std::convertible_to<double>;
  { m.m22 } -> std::convertible_to<double>;
  { m.m23 } -> std::convertible_to<double>;
  { m.m31 } -> std::convertible_to<double>;
  { m.m32 } -> std::convertible_to<double>;
  { m.m33 } -> std::convertible_to<double>;
};

struct Vec3dView {
  double& x;
  double& y;
  double& z;
};

template <size_t I, Vec3dLike Vec> inline constexpr double& at(Vec& v) {
  static_assert(I >= 0 && I < 3, "Invalid indexing of Vec3d");
  if constexpr (I == 0)
    return &v.x;
  else if constexpr (I == 1)
    return &v.y;
  else if constexpr (I == 2)
    return &v.z;
  else
    PANIC("");
}

constexpr inline double& at(Vec3d& v, size_t index) {
  double* base = &v.x;
  return base[index];
}

template <size_t I, size_t J, Mat3dLike Matrix> constexpr inline double& at(Matrix& m) {
  static_assert(I >= 0 && I < 3);
  static_assert(J >= 0 && J < 3);
  if constexpr (I == 0 && J == 0)
    return m.m11;
  else if constexpr (I == 0 && J == 1)
    return m.m12;
  else if constexpr (I == 0 && J == 2)
    return m.m13;
  else if constexpr (I == 1 && J == 0)
    return m.m21;
  else if constexpr (I == 1 && J == 1)
    return m.m22;
  else if constexpr (I == 1 && J == 2)
    return m.m23;
  else if constexpr (I == 2 && J == 0)
    return m.m31;
  else if constexpr (I == 2 && J == 1)
    return m.m32;
  else if constexpr (I == 2 && J == 2)
    return m.m33;
  else
    static_assert(std::false_type{});
}

template <size_t I> inline constexpr Vec3dView column(Mat3d& m) {
  static_assert(I >= 0 && I < 3, "Invalid column indexing.");
  if constexpr (I == 0)
    return Vec3dView{m.m11, m.m21, m.m31};
  if constexpr (I == 1)
    return Vec3dView{m.m12, m.m22, m.m32};
  if constexpr (I == 2)
    return Vec3dView{m.m13, m.m23, m.m33};
};

inline Vec3dView column(size_t i, Mat3d& m) {
  double* ptr = &m.m11;
  return Vec3dView{ptr[i], ptr[i + 3], ptr[i + 6]};
}

inline Vec3d column(size_t i, const Mat3d& m) {
  const double* ptr = &m.m11;
  return Vec3d{ptr[i], ptr[i + 3], ptr[i + 6]};
}

template <size_t I> inline constexpr Vec3dView row(Mat3d& m) {
  static_assert(I >= 0 && I < 3, "Invalid raw indexing.");
  if constexpr (I == 0)
    return Vec3dView{m.m11, m.m12, m.m13};
  if constexpr (I == 1)
    return Vec3dView{m.m21, m.m22, m.m23};
  if constexpr (I == 2)
    return Vec3dView{m.m31, m.m32, m.m33};
};

inline Vec3dView row(size_t i, Mat3d& m) {
  double* ptr = &m.m11;
  size_t index = 3 * i;
  return Vec3dView{ptr[index], ptr[index + 1], ptr[index + 2]};
}

inline Vec3d row(size_t i, const Mat3d& m) {
  const double* ptr = &m.m11;
  size_t index = 3 * i;
  return Vec3d{ptr[index], ptr[index + 1], ptr[index + 2]};
}

template <Vec3dLike Vec> inline void scale(Vec& v, double scale) {
  v.x *= scale;
  v.y *= scale;
  v.z *= scale;
}

} // namespace matrix_traits

namespace mem_size {

constexpr size_t KB = 1 << 10;
constexpr size_t MB = 1 << 20;
constexpr size_t GB = 1 << 30;

constexpr size_t READ_BUFFER_SIZE = 1 * MB;
constexpr size_t FILE_BUFFER_SIZE = READ_BUFFER_SIZE;
constexpr size_t LZMA_STREAM_BUF_SIZE = READ_BUFFER_SIZE;
constexpr size_t LZMA_INTERNAL_BUF_SIZE = READ_BUFFER_SIZE;
constexpr size_t BZ2_STREAM_BUF_SIZE = READ_BUFFER_SIZE;
constexpr size_t BZ2_INTERNAL_BUF_SIZE = READ_BUFFER_SIZE;

constexpr size_t WRITE_BUF_SIZE = 1 * MB;
// constexpr size_t WRITE_BUF_THRESH = FILE_BUFFER_SIZE - 100 * KB;
constexpr size_t WRITE_BUF_THRESH = 1 * KB;
} // namespace mem_size

namespace write_utils {

template <typename S>
concept WriteSink = requires(S& s, const char* p, size_t n) {
  { s.write(p, n) } -> std::same_as<void>;
};

struct FloatFormat {
  std::chars_format fmt = std::chars_format::fixed;
  int precision = 6;
  int width = 0;
  char fill = ' ';
  bool align_right = true;
};

template <size_t BUF_SIZE, WriteSink Sink>
  requires(BUF_SIZE > mem_size::WRITE_BUF_THRESH)
class IOBuffer {
  static constexpr size_t TOCHARS_BUF = 64;
  static constexpr size_t FLUSH_TRESHOLD = BUF_SIZE - mem_size::WRITE_BUF_THRESH;

  std::array<char, BUF_SIZE> m_buf{};
  size_t m_pos = 0;
  Sink& m_sink; // abstract sink

public:
  FloatFormat format{};

  explicit IOBuffer(Sink& sink) : m_sink(sink) {}
  ~IOBuffer() { flush(); }

  IOBuffer(const IOBuffer&) = delete;
  IOBuffer& operator=(const IOBuffer&) = delete;
  IOBuffer(IOBuffer&&) = delete;
  IOBuffer& operator=(IOBuffer&&) = delete;

  inline size_t size() const { return m_pos; }
  inline const char* data() const { return m_buf.data(); }

  void clear() { m_pos = 0; }

  void flush() {
    if (m_pos > 0) {
      m_sink.write(m_buf.data(), m_pos);
      m_pos = 0;
    }
  }

  // reserver space, flush if needed
  inline void ensure_space(size_t n) {
    if (cmp::gt(m_pos + n, FLUSH_TRESHOLD)) [[unlikely]]
      flush();
  }

  inline void insert_char(char value) {
    ensure_space(1);
    m_buf[m_pos++] = value;
  }

  inline void insert_string(std::string_view value) {
    ensure_space(value.size());
    std::memcpy(m_buf.data() + m_pos, value.data(), value.size());
    m_pos += value.size();
  }

  template <std::floating_point T> inline void insert_float(T value, const FloatFormat& fmt) {
    char tmp[TOCHARS_BUF];
    auto [ptr, ec] = std::to_chars(tmp, tmp + TOCHARS_BUF, value, fmt.fmt, fmt.precision);
    if (ec != std::errc()) [[unlikely]]
      return;
    insert_padded(tmp, ptr - tmp, fmt);
  }

  template <std::integral T> inline void insert_int(T value, const FloatFormat& fmt) {
    char tmp[TOCHARS_BUF];
    auto [ptr, ec] = std::to_chars(tmp, tmp + TOCHARS_BUF, value);
    if (ec != std::errc()) [[unlikely]]
      return;
    insert_padded(tmp, ptr - tmp, fmt);
  }

  template <matrix_traits::Vec3dLike V> inline void insert_vec3d(const V& v, const FloatFormat& fmt) {
    insert_float(v.x, fmt);
    space();
    insert_float(v.y, fmt);
    space();
    insert_float(v.z, fmt);
  }

  inline void insert_padded(const char* src, size_t len, const FloatFormat& fmt) {
    size_t pad = (fmt.width > 0 && len < static_cast<size_t>(fmt.width)) ? fmt.width - len : 0;
    size_t total = pad + len;
    ensure_space(total);

    if (fmt.align_right && pad) {
      std::fill_n(m_buf.data() + m_pos, pad, fmt.fill);
      std::memcpy(m_buf.data() + m_pos + pad, src, len);
    } else {
      std::memcpy(m_buf.data() + m_pos, src, len);
      if (pad)
        std::fill_n(m_buf.data() + m_pos + len, pad, fmt.fill);
    }
    m_pos += total;
  }

  inline void space() { insert_char(' '); }
  inline void newline() { insert_char('\n'); }

  // generic dispatcher
  template <typename T> IOBuffer& insert(const T& value) {
    if constexpr (std::is_same_v<T, char>) {
      insert_char(value);
    } else if constexpr (std::is_same_v<T, std::string_view>) {
      insert_string(value);
    } else if constexpr (matrix_traits::Vec3dLike<T>) {
      insert_vec3d(value, format);
    } else if constexpr (std::is_floating_point_v<T>) {
      insert_float(value, format);
    } else if constexpr (std::is_integral_v<T>) {
      insert_int(value, format);
    } else if constexpr (std::is_convertible_v<T, std::string_view>) {
      insert_string(std::string_view(value));
    } else {
      static_assert([] { return false; }(), "Unsupported type");
    }
    return *this;
  }

  template <typename T> void insert(const T& value, const FloatFormat& fmt) {
    if constexpr (matrix_traits::Vec3dLike<T>) {
      insert_vec3d(value, fmt);
    } else if constexpr (std::is_floating_point_v<T>) {
      insert_float(value, fmt);
    } else if constexpr (std::is_integral_v<T>) {
      insert_int(value, fmt);
    } else {
      static_assert([] { return false; }(), "Unsupported type");
    }
  }
};
} // namespace write_utils

namespace io {

using namespace tokens;

template <policy::ValidPolicy P> struct IOContextTmpl {

  using DataBuffer = typename P::ContextDataBuffer;
  using SpeciesMap = typename base_types::SpeciesMap;

  enum ContextFlags : uint16_t {
    REMAP_ATOM = 1 << 1,
    DOMAIN_ONLY = 1 << 2,
    TRICLINIC = 1 << 3,
    NO_VELOCITY = 1 << 4,
  };

  template <uint16_t Flag> inline bool has() { return (flags & (Flag)); }

  template <int16_t Flag> void set(bool toggle) { (toggle) ? flags |= Flag : flags &= Flag; }

  uint16_t flags = 0;
  P::Vec3d timer{};
  DataBuffer data{};
  SpeciesMap species{};
};

using IOContext = IOContextTmpl<policy::Traits>;

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
    return "";
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
  TextFileHandler(std::string filepath) : m_path(std::move(filepath)) {}
  virtual ~TextFileHandler() = default;
  virtual void clear() noexcept = 0;
  virtual void seek(uint64_t position) = 0;
  virtual uint64_t tell(void) = 0;
  virtual size_t read(char* buffer, size_t size) = 0;
  virtual size_t write(const char* data, size_t size) = 0;

protected:
  const std::string& path() { return m_path; };
  std::string m_path;
};

class ASCIIFileHandler final : public TextFileHandler {
public:
  ASCIIFileHandler(std::string filepath, FileMode mode) : TextFileHandler(std::move(filepath)) {
    m_fptr = std::fopen(path().c_str(), convert_mode_to_char(mode));
    if (!m_fptr)
      PANIC("Could not open the file: %s", path().c_str());
  }

  ~ASCIIFileHandler() override {
    if (m_fptr)
      std::fclose(m_fptr);
  }

  inline void clear() noexcept override { std::clearerr(m_fptr); };

  void seek(uint64_t pos) override {
    static_assert(sizeof(uint64_t) == sizeof(off64_t));
    off64_t ierr = fseeko(m_fptr, static_cast<off64_t>(pos), SEEK_SET);
    if (ierr != 0)
      PANIC("Fail to seek inside file.");
  }

  uint64_t tell() override { return static_cast<uint64_t>(ftello(m_fptr)); }

  size_t read(char* data, size_t size) noexcept override {
    size_t n = std::fread(data, 1, size, m_fptr);

    if (n < size && std::feof(m_fptr))
      return n;

    if (std::ferror(m_fptr))
      PANIC("Something when reading the file %s", path().c_str());

    return n;
  }

  size_t write(const char* data, size_t size) noexcept override {
    size_t n = std::fwrite(data, 1, size, m_fptr);
    if (std::ferror(m_fptr) != 0)
      PANIC("Fail to write file %s", path().c_str());
    return n;
  }

private:
  std::FILE* m_fptr = nullptr;
};

#ifdef USE_ZLIB

class GzipFileHandler final : public TextFileHandler {
public:
  GzipFileHandler(std::string path, FileMode mode) : TextFileHandler(std::move(path)) {
    m_fptr = gzopen64(m_path.c_str(), mode_to_char(mode));
    if (!m_fptr)
      PANIC("Unable to open gzip file: %s", m_path.c_str());
  }

  ~GzipFileHandler() override {
    if (m_fptr)
      gzclose(m_fptr);
  }

  void clear() noexcept override { gzclearerr(m_fptr); }

  void seek(uint64_t cursor) override {
    static_assert(sizeof(uint64_t) == sizeof(z_off64_t));
    if (gzseek64(m_fptr, static_cast<z_off64_t>(cursor), SEEK_SET) == -1) {
      PANIC("Error while decompressing gzip file : %s", last_gz_error());
    }
  }

  uint64_t tell() override {
    static_assert(sizeof(uint64_t) == sizeof(z_off64_t));
    z_off64_t pos = gztell64(m_fptr);
    if (pos < 0)
      PANIC("Error telling position in gzip file: %s", last_gz_error());
    return static_cast<uint64_t>(pos);
  }

  size_t read(char* buffer, size_t size) override {
    int n = gzread(m_fptr, buffer, safe_cast(size));
    if (n < 0)
      PANIC("Error while reading gzip file : %s", last_gz_error());
    return static_cast<size_t>(n);
  }

  size_t write(const char* buffer, size_t size) override {
    int n = gzwrite(m_fptr, buffer, safe_cast(size));
    if (n == 0)
      PANIC("Error writing gzip file: %s", last_gz_error());
    return static_cast<size_t>(n);
  }

private:
  const char* mode_to_char(FileMode mode) {
    switch (mode) {
    case FileMode::READ:
      return "rb";
    case FileMode::APPEND:
      return "ab"; // supported by zlib
    case FileMode::WRITE:
      return "wb";
    default:
      PANIC("Unsupported FileMode for GzipFileHandler");
    }
  }

  const char* last_gz_error() const noexcept {
    int status = Z_OK;
    const char* msg = gzerror(m_fptr, &status);
    return (status == Z_OK ? "Unknow ZLIB Error" : msg);
  }

  static unsigned safe_cast(size_t value) {
    constexpr size_t max = std::numeric_limits<unsigned>::max();
    if (cmp::gt(value, max))
      PANIC("%zu is too big for zlib API (max = %u)", value, max);
    return static_cast<unsigned>(value);
  }

  gzFile m_fptr = nullptr;
};

#endif

#ifdef USE_BZIP2

class Bzip2FileHandler final : public TextFileHandler {
public:
  Bzip2FileHandler(std::string filepath, FileMode mode) : TextFileHandler(std::move(filepath)) {

    const char* m = nullptr;
    switch (mode) {
    case FileMode::READ:
      m = "rb";
      m_end_bz2_stream = &BZ2_bzDecompressEnd;
      check_bz2_retcode(BZ2_bzDecompressInit(&m_bz2_stream, 0, 0));
      m_read_mode = true;
      break;

    case FileMode::WRITE:
      m = "wb";
      m_end_bz2_stream = &BZ2_bzCompressEnd;
      check_bz2_retcode(BZ2_bzCompressInit(&m_bz2_stream, 9, 0, 30));
      m_read_mode = true;
      break;

    default:
      PANIC("Unsupported mode for Bzip2 format");
      break;
    }

    m_fptr = std::fopen(path().c_str(), m);
    if (!m_fptr) {
      if (m_end_bz2_stream)
        m_end_bz2_stream(&m_bz2_stream);
      PANIC("Unable to open bz2 file: %s", m_path.c_str());
    }
  }

  ~Bzip2FileHandler() override {
    // flush pending compressed data on write
    if (!m_read_mode) {
      int status;
      do {
        m_bz2_stream.next_out = m_bz2_buffer.data();
        m_bz2_stream.avail_out = m_bz2_buffer.size();

        status = BZ2_bzCompress(&m_bz2_stream, BZ_FINISH);
        check_bz2_retcode(status);

        size_t remain = m_bz2_buffer.size() - m_bz2_stream.avail_out;
        if (remain > 0 && std::fwrite(m_bz2_buffer.data(), 1, remain, m_fptr) != remain)
          PANIC("bzip2: failed to write compressed tail");
      } while (status != BZ_STREAM_END);
    }

    // close stream / file
    if (m_end_bz2_stream)
      m_end_bz2_stream(&m_bz2_stream);
    if (m_fptr)
      std::fclose(m_fptr);
  }

  void clear() noexcept override { std::clearerr(m_fptr); }
  uint64_t tell() override { return m_uncompressed_pos; }

  void seek(uint64_t cursor) override {

    // bzip2 is a stream based compression format, so random access to file is not supported.
    // inefficient implementation: re-decompressing the file from the begining
    // not really an issue if we read from the begining of the file.

    if (!m_read_mode) {
      PANIC("seek is not supported in write mode for bzip2");
    }

    if (m_end_bz2_stream)
      m_end_bz2_stream(&m_bz2_stream);
    check_bz2_retcode(BZ2_bzDecompressInit(&m_bz2_stream, 0, 0));
    std::fseek(m_fptr, 0, SEEK_SET);
    m_uncompressed_pos = 0;

    char tmpbuf[mem_size::BZ2_STREAM_BUF_SIZE];
    while (cursor > 0) {
      size_t toread = std::min<size_t>(mem_size::BZ2_STREAM_BUF_SIZE, cursor);
      size_t got = this->read(tmpbuf, toread);
      if (got == 0)
        break; // EOF
      cursor -= got;
    }
  }

  size_t read(char* buffer, size_t size) override {
    if (!m_read_mode) {
      PANIC("Attempt to read a bzip2 file open in write mode.");
    }

    m_bz2_stream.next_out = buffer;
    m_bz2_stream.avail_out = safe_cast(size);

    while (m_bz2_stream.avail_out > 0) {
      if (m_bz2_stream.avail_in == 0) {
        size_t n = std::fread(m_bz2_buffer.data(), 1, m_bz2_buffer.size(), m_fptr);
        if (n == 0) {
          if (std::feof(m_fptr))
            break;
          PANIC("Error while reading bzip2 file: %s", path().c_str());
        }
        m_bz2_stream.next_in = m_bz2_buffer.data();
        m_bz2_stream.avail_in = safe_cast(n);
      }

      int code = BZ2_bzDecompress(&m_bz2_stream);
      if (code == BZ_STREAM_END)
        break;
      check_bz2_retcode(code);
    }

    size_t nread = static_cast<size_t>(size - m_bz2_stream.avail_out);
    m_uncompressed_pos += nread;
    return nread;
  }

  size_t write(const char* data, size_t size) override {
    if (m_read_mode) {
      PANIC("Attempt to write to bzip2 file in read mode");
    }

    m_bz2_stream.next_in = const_cast<char*>(data);
    m_bz2_stream.avail_in = safe_cast(size);

    while (m_bz2_stream.avail_in > 0) {
      m_bz2_stream.next_out = m_bz2_buffer.data();
      m_bz2_stream.avail_out = m_bz2_buffer.size();

      check_bz2_retcode(BZ2_bzCompress(&m_bz2_stream, BZ_RUN));

      size_t remain = m_bz2_buffer.size() - m_bz2_stream.avail_out;
      if (remain > 0 && std::fwrite(m_bz2_buffer.data(), 1, remain, m_fptr) != remain) {
        PANIC("Error writing bzip2 compressed data to file: %s", path().c_str());
      }
    }
    m_uncompressed_pos += size;
    return size;
  }

  unsigned safe_cast(uint64_t size) {
    constexpr uint64_t max = std::numeric_limits<unsigned>::max();
    if (size > max) {
      PANIC("%llu is too large for bzlib unsigned parameter", (unsigned long long)size);
    }
    return static_cast<unsigned>(size);
  }

  inline void check_bz2_retcode(int code) {
    switch (code) {
    case BZ_OK:
    case BZ_RUN_OK:
    case BZ_FLUSH_OK:
    case BZ_FINISH_OK:
    case BZ_STREAM_END:
      return;
    default:
      PANIC("bzip2 error (code=%d)", code);
    }
  }

private:
  std::FILE* m_fptr = nullptr;
  bz_stream m_bz2_stream{};
  int (*m_end_bz2_stream)(bz_stream*) = nullptr;
  bool m_read_mode{};
  std::array<char, mem_size::BZ2_INTERNAL_BUF_SIZE> m_bz2_buffer;
  uint64_t m_uncompressed_pos = 0;
};

#endif

#ifdef USE_LZMA

class XzFileHandler final : public TextFileHandler {
public:
  XzFileHandler(std::string filepath, FileMode mode) : TextFileHandler(std::move(filepath)) {

    const char* m = nullptr;
    switch (mode) {
    case FileMode::READ:
      m = "rb";
      m_xz_stream = LZMA_STREAM_INIT;
      start_lzma_decoder_stream(&m_xz_stream);
      m_read_mode = true;
      break;

    case FileMode::WRITE:
      m = "wb";
      m_xz_stream = LZMA_STREAM_INIT;
      start_lzma_encoder_stream(&m_xz_stream);
      m_read_mode = false;
      break;

    default:
      PANIC("Unsupported mode for XZ format");
      break;
    }

    m_fptr = std::fopen(path().c_str(), m);
    if (!m_fptr) {
      lzma_end(&m_xz_stream);
      PANIC("Unable to open xz file: %s", m_path.c_str());
    }
  }

  ~XzFileHandler() override {
    if (!m_read_mode)
      flush_encoder();

    lzma_end(&m_xz_stream);
    if (m_fptr)
      std::fclose(m_fptr);
  }

  void clear() noexcept override { std::clearerr(m_fptr); }
  uint64_t tell() override { return m_uncompressed_pos; }

  void seek(uint64_t cursor) override {
    if (!m_read_mode) {
      PANIC("Seek not supported in write mode for xz files.");
    }
    // Re-decompress from beginning
    lzma_end(&m_xz_stream);
    m_xz_stream = LZMA_STREAM_INIT;
    start_lzma_decoder_stream(&m_xz_stream);

    std::fseek(m_fptr, 0, SEEK_SET);
    char buf[mem_size::LZMA_INTERNAL_BUF_SIZE];

    while (cursor > sizeof(buf)) {
      size_t n = read(buf, sizeof(buf));
      assert(n == sizeof(buf));
      cursor -= n;
    }

    [[maybe_unused]] size_t n = read(buf, static_cast<size_t>(cursor));
    assert(n == cursor);
    m_uncompressed_pos = cursor;
  }

  size_t read(char* buffer, size_t size) override {

    if (!m_read_mode)
      PANIC("Attempt to read an xz file opened in write mode.");

    m_xz_stream.next_out = reinterpret_cast<uint8_t*>(buffer);
    m_xz_stream.avail_out = size;

    while (m_xz_stream.avail_out > 0) {
      if (m_xz_stream.avail_in == 0) {
        size_t n = std::fread(m_xz_buffer.data(), 1, m_xz_buffer.size(), m_fptr);
        if (n == 0) {
          if (std::feof(m_fptr))
            break;
          PANIC("Error while reading xz file: %s", m_path.c_str());
        }
        m_xz_stream.next_in = m_xz_buffer.data();
        m_xz_stream.avail_in = n;
      }

      lzma_action action = (std::feof(m_fptr) && m_xz_stream.avail_in == 0) ? LZMA_FINISH : LZMA_RUN;
      lzma_ret ret = lzma_code(&m_xz_stream, action);

      if (ret == LZMA_STREAM_END)
        break;

      // Output buffer full, not an actual error
      if (ret == LZMA_BUF_ERROR && m_xz_stream.avail_out == 0)
        break;

      check_lzma_ret(ret);
    }

    size_t nread = size - m_xz_stream.avail_out;
    m_uncompressed_pos += nread;
    return nread;
  }

  size_t write(const char* data, size_t size) override {

    if (m_read_mode)
      PANIC("Attempt to write to an xz file opened in read mode.");

    m_xz_stream.next_in = reinterpret_cast<const uint8_t*>(data);
    m_xz_stream.avail_in = size;

    while (m_xz_stream.avail_in > 0) {
      m_xz_stream.next_out = m_xz_buffer.data();
      m_xz_stream.avail_out = m_xz_buffer.size();

      lzma_ret ret = lzma_code(&m_xz_stream, LZMA_RUN);
      if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
        check_lzma_ret(ret);
      }

      size_t n = m_xz_buffer.size() - m_xz_stream.avail_out;
      if (n > 0 && std::fwrite(m_xz_buffer.data(), 1, n, m_fptr) != n) {
        PANIC("Error writing compressed data to file: %s", m_path.c_str());
      }
    }

    m_uncompressed_pos += size;
    return size;
  }

  static void check_lzma_ret(lzma_ret ret) {
    switch (ret) {
    case LZMA_OK:
    case LZMA_STREAM_END:
      return;
    case LZMA_NO_CHECK:
    case LZMA_GET_CHECK:
      PANIC("lzma: no integrity check.");
    case LZMA_MEM_ERROR:
    case LZMA_MEMLIMIT_ERROR:
      PANIC("lzma: memory allocation failed.");
    case LZMA_FORMAT_ERROR:
      PANIC("lzma: invalid .xz format.");
    case LZMA_OPTIONS_ERROR:
      PANIC("lzma: invalid compression options.");
    case LZMA_DATA_ERROR:
      PANIC("lzma: corrupted or truncated file.");
    case LZMA_UNSUPPORTED_CHECK:
      PANIC("lzma: unsupported integrity check.");
    case LZMA_PROG_ERROR:
      PANIC("lzma: library internal bug.");
    case LZMA_BUF_ERROR:
      // handled in read/write loop
      return;
    default:
      PANIC("Uknown lzma error code :/");
    }
  }

  void flush_encoder() {
    if (m_read_mode)
      return;

    m_xz_stream.next_in = nullptr;
    m_xz_stream.avail_in = 0;

    while (true) {
      m_xz_stream.next_out = m_xz_buffer.data();
      m_xz_stream.avail_out = m_xz_buffer.size();

      lzma_ret ret = lzma_code(&m_xz_stream, LZMA_FINISH);

      size_t n = m_xz_buffer.size() - m_xz_stream.avail_out;
      if (n > 0 && std::fwrite(m_xz_buffer.data(), 1, n, m_fptr) != n) {
        PANIC("Error flushing compressed data to file: %s", m_path.c_str());
      }

      if (ret == LZMA_STREAM_END)
        break;
      check_lzma_ret(ret);
    }
  }

  static void start_lzma_decoder_stream(lzma_stream* stream) {
    constexpr uint64_t memlimit = std::numeric_limits<uint64_t>::max();
    check_lzma_ret(lzma_stream_decoder(stream, memlimit, LZMA_TELL_UNSUPPORTED_CHECK | LZMA_CONCATENATED));
  }

  static void start_lzma_encoder_stream(lzma_stream* stream) {
    uint32_t preset = 6 | LZMA_PRESET_EXTREME; // compression level
    check_lzma_ret(lzma_easy_encoder(stream, preset, LZMA_CHECK_CRC64));
  }

private:
  std::FILE* m_fptr = nullptr;
  lzma_stream m_xz_stream{};
  std::array<uint8_t, mem_size::LZMA_INTERNAL_BUF_SIZE> m_xz_buffer;
  uint64_t m_uncompressed_pos = 0;
  bool m_read_mode{};
};

#endif

namespace {

static std::unique_ptr<TextFileHandler> create_file_handler(const std::string& path, FileMode mode,
                                                            FileCompression compression) {
  switch (compression) {
  case FileCompression::NONE:
    return std::make_unique<ASCIIFileHandler>(path, mode);
    break;
  case FileCompression::GZIP:
    if constexpr (has_zlib)
      return std::make_unique<GzipFileHandler>(path, mode);
    break;
  case FileCompression::BZIP2:
    if constexpr (has_bzip2)
      return std::make_unique<Bzip2FileHandler>(path, mode);
    break;
  case FileCompression::XZ:
    if constexpr (has_lzma)
      return std::make_unique<XzFileHandler>(path, mode);
    break;
  default:
    break;
  }
  PANIC("Unable to instantiate a file handler.");
  return nullptr;
}

} // namespace

/* ------------------------------------------------------------------------- */

class File {
protected:
  std::string m_path;
  FileMode m_mode;
  FileCompression m_compression;
  bool m_eof = false;

  File(const std::string& path, FileMode mode, FileCompression compression)
      : m_path(std::move(path)), m_mode(mode), m_compression(compression) {}

public:
  inline const std::string& path() const noexcept { return m_path; };
  inline FileMode mode() const noexcept { return m_mode; };
  inline FileCompression compression() const noexcept { return m_compression; };
  inline bool eof() const noexcept { return m_eof; }

  virtual std::string_view get_line() = 0;
  virtual void seek(uint64_t) = 0;
  virtual uint64_t tell() = 0;
  virtual void reset() = 0;
};

class MemoryMappedFile : public File {
public:
  MemoryMappedFile(const std::string& filepath) : File(std::move(filepath), FileMode::READ, FileCompression::NONE) {

    m_fd = open(path().c_str(), O_RDONLY);

    if (m_fd < 0)
      PANIC("Unable to open file: %s", path().c_str());

    struct stat fst{};
    if (fstat(m_fd, &fst) != 0)
      PANIC("Unable to fetch file stats for %s", path().c_str());

    m_size = fst.st_size;
    m_data = static_cast<const char*>(mmap(nullptr, m_size, PROT_READ, MAP_PRIVATE, m_fd, 0));
    if (m_data == nullptr)
      PANIC("Unable to memory map file: %s", path().c_str());

    m_ptr = m_data;
    m_end = m_data + m_size;
  }

  ~MemoryMappedFile() {
    munmap(const_cast<char*>(m_data), m_size);
    close(m_fd);
  }

  inline void seek(uint64_t pos) override {
    if (pos >= m_size)
      PANIC("Invalid seek position.");
    m_eof = false;
    m_ptr = static_cast<const char*>(m_data + pos);
  }

  inline uint64_t tell() override { return static_cast<uint64_t>(m_ptr - m_data); }

  inline void reset() override {
    m_eof = false;
    seek(0);
  }

  std::string_view get_line() override {

    if (m_ptr >= m_end) {
      m_eof = true;
      return {};
    }

    const char* start = m_ptr;
    const char* new_line = static_cast<const char*>(memchr(start, '\n', m_end - start));

    // 0xFFFFFFFF... if found, 0 otherwise
    const uintptr_t mask = -(new_line != nullptr);
    const char* end = reinterpret_cast<const char*>((reinterpret_cast<uintptr_t>(m_end) & ~mask) |
                                                    (reinterpret_cast<uintptr_t>(new_line) & mask));
    size_t len = static_cast<size_t>(end - start);
    m_ptr = end + (new_line != nullptr);

    return {start, len};
  }

private:
  int m_fd;
  size_t m_size;
  const char* m_data;
  const char* m_ptr;
  const char* m_end;
};

class TextFile : public File {
private:
  std::array<char, mem_size::FILE_BUFFER_SIZE> m_buf{};
  const char* m_line_start; // current read cursor in buffer
  const char* m_buf_end;    // one-past-last valid byte in buffer

  uint64_t m_pos = 0;      // absolute file position of buffer start
  bool m_hdlr_eof = false; // whether handler reported EOF
  std::unique_ptr<TextFileHandler> m_handler;

  inline bool has_buffer_data() const noexcept { return m_line_start < m_buf_end; }
  inline bool is_buffer_initialized() const noexcept { return m_buf_end != m_buf.data(); }

  void fill_buffer(size_t offset = 0) {

    if (offset >= m_buf.size()) {
      PANIC("Offset value larger than buffer size (%ld >> %ld)!!!", offset, m_buf.size());
    }

    size_t nread = m_buf.size() - offset;
    size_t n = m_handler->read(m_buf.data() + offset, nread);
    m_buf_end = m_buf.data() + offset + n;

    // if the number of bytes read is less that the asked value we reached the EOF
    if (n < nread)
      m_hdlr_eof = true;

    if (offset == 0) {
      m_pos = m_handler->tell() - n;
    }

    // we advanced in the file this many bytes
    // m_pos += n;
  };

public:
  TextFile(const std::string& filepath, FileMode mode, FileCompression compression)
      : File(std::move(filepath), mode, compression), m_line_start(m_buf.data()), m_buf_end(m_buf.data()),
        m_handler(create_file_handler(filepath, mode, compression)) {}

  uint64_t tell() override {
    uint64_t delta = static_cast<uint64_t>(m_line_start - m_buf.data());
    return m_pos + delta;
  };

  void seek(uint64_t pos) override {
    m_handler->seek(pos);
    m_pos = pos;
    m_line_start = m_buf.data();
    m_buf_end = m_buf.data();
    m_hdlr_eof = false;
    m_eof = false;
  }

  void clear() noexcept {
    m_eof = false;
    m_hdlr_eof = false;
    m_line_start = m_buf.data();
    m_buf_end = m_buf.data();
    m_handler->clear();
  }

  void reset() override {
    clear();
    seek(0);
  }

  std::string_view get_line() override {

    if (m_eof)
      return {};

    if (!is_buffer_initialized()) {
      fill_buffer(0);
      m_line_start = m_buf.data();

      if (!has_buffer_data() && m_hdlr_eof) {
        m_eof = true;
        return {};
      }
    }

    while (true) {

      if (!has_buffer_data()) {
        fill_buffer(0);
        m_line_start = m_buf.data();

        if (!has_buffer_data()) {
          m_eof = true;
          return {};
        }
      }

      size_t avail = static_cast<size_t>(m_buf_end - m_line_start);
      const char* newline = static_cast<const char*>(std::memchr(m_line_start, '\n', avail));

      if (newline) {
        size_t len = static_cast<size_t>(newline - m_line_start);
        std::string_view line{m_line_start, len};
        m_line_start = static_cast<const char*>(newline + 1);
        return line;
      }

      if (m_hdlr_eof) {
        m_eof = true;
        std::string_view line{m_line_start, avail};
        m_line_start = m_buf_end;
        return line;
      }

      // no newline, not EOF yes -> shift remainder and refill
      if (avail > 0) {
        std::memmove(m_buf.data(), m_line_start, avail);
        m_line_start = m_buf.data();
        fill_buffer(avail);
      }
    }
  }

  inline void write(const char* data, size_t n) { m_handler->write(data, n); }
};

/* ------------------------------------------------------------------------- */

struct FormatMetadata {
  const std::string name = "";
  const std::string description = "";
  const std::vector<std::string> aliases{};
  size_t alias_count = 0;
};

template <typename Derived> struct IOTraitsBase {
  static FormatMetadata metadata() {
    return {
        .name = std::string(Derived::name),
        .description = std::string(Derived::description),
        .aliases = std::vector<std::string>(Derived::aliases.begin(), Derived::aliases.end()),
    };
  }
};

template <typename T>
concept HasValidMetadataTraits = requires {
  typename T::Traits;

  { T::Traits::name } -> std::convertible_to<const char*>;
  {
    T::Traits::aliases
  } -> std::convertible_to<const std::array<const char*, std::tuple_size<decltype(T::Traits::aliases)>::value>&>;
  { T::Traits::description } -> std::convertible_to<const char*>;
  { T::Traits::metadata() } -> std::convertible_to<FormatMetadata>;

  requires T::Traits::name[0] != '\0';
  requires T::Traits::description[0] != '\0';
  requires T::Traits::aliases[0][0] != '\0';
};

template <typename Format>
  requires HasValidMetadataTraits<Format>
const FormatMetadata GetFormatMetadata() {
  return Format::Traits::metadata();
}

/* ------------------------------------------------------------------------- */

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

/* ------------------------------------------------------------------------- */

class TextParser : public Parser {

protected:
  bool m_eof = false;
  std::vector<size_t> m_loc{};
  // std::unique_ptr<File> m_file = nullptr;

  TextFile m_file;

  std::array<std::string_view, 10> m_line_buffer{};
  TokenSet m_tokens{};

  inline std::string_view& current_line(void) { return m_line_buffer[0]; };
  inline const std::string_view& current_line(void) const { return m_line_buffer[0]; }

public:
  TextParser(std::string filepath, FileMode mode, FileCompression compression)
      : m_file(std::move(filepath), mode, compression) {
    // if (!details::file_exists(filepath.c_str()))
    //   PANIC("File do not exists !!!");
    // m_file = std::make_unique<TextFile>(std::move(filepath), mode, compression);
    // m_file = std::make_unique<MemoryMappedFile>(filepath);
    // uint64_t fsize = details::file_size(filepath.c_str());
    // if (fsize < 64 * mem_size::MB && mode == FileMode::READ &&
    //     compression == FileCompression::NONE) {
    //   m_file = std::make_unique<MemoryMappedFile>(filepath);
    // } else {
    //   m_file = std::make_unique<TextFile>(std::move(filepath), mode, compression);
    // }
  }

  virtual ~TextParser() override = default;

  virtual int64_t next() = 0;
  virtual bool parse(IOContext& ctx) = 0;

  inline bool eof() { return m_eof; }
  // TODO:: void scan();

  bool operator()(IOContext& ctx) override {
    std::chrono::time_point tstart = std::chrono::high_resolution_clock::now();

    m_file.reset();
    bool success = this->parse(ctx);

    std::chrono::time_point tend = std::chrono::high_resolution_clock::now();
    ctx.timer = {
        static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count()),
        static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(tend - tstart).count()),
        static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(tend - tstart).count()),
    };

    return success;
  }

  // TODO:
  bool operator()(IOContext& ctx, size_t /*index*/) override { return (*this)(ctx); };
  inline size_t size() override { return m_loc.size(); };
};

/* ------------------------------------------------------------------------- */

class Writer {
public:
  Writer() = default;
  virtual ~Writer() = default;
  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;
  Writer(Writer&&) = delete;
  Writer& operator=(Writer&&) = delete;

  virtual bool write(IOContext& ctx) = 0;
  virtual bool operator()(io::IOContext&) { return false; };
};

/* ------------------------------------------------------------------------- */

class TextWriter : public Writer {
protected:
  TextFile m_file;
  write_utils::IOBuffer<mem_size::WRITE_BUF_SIZE, TextFile> m_buf;

public:
  TextWriter(std::string filepath, FileMode mode, FileCompression compression)
      : m_file(std::move(filepath), mode, compression), m_buf(m_file) {}

  virtual ~TextWriter() override = default;

  bool operator()(IOContext& ctx) override {
    std::chrono::time_point tstart = std::chrono::high_resolution_clock::now();

    m_file.reset();
    bool success = this->write(ctx);

    std::chrono::time_point tend = std::chrono::high_resolution_clock::now();
    ctx.timer = {
        static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count()),
        static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(tend - tstart).count()),
        static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(tend - tstart).count()),
    };

    return success;
  }
};

/* ------------------------------------------------------------------------- */

namespace factory {

template <typename I> struct IOFactory {
  using creator_t = std::unique_ptr<I> (*)(std::string, FileMode, FileCompression);

  struct Format {
    const FormatMetadata metadata = FormatMetadata{};
    const creator_t creator = nullptr;
    bool is_valid = true;
  };

  std::unordered_map<std::string, size_t> alias_to_index;
  std::unordered_map<size_t, Format> index_to_parser;
  std::array<Format, 32> formats;
  size_t size = 0;

  IOFactory() {}

  template <typename T>
  static std::unique_ptr<I> create_instance(std::string filepath, FileMode mode, FileCompression compression) {
    return std::make_unique<T>(std::move(filepath), mode, compression);
  }

  template <typename T> void register_format() {
    size_t index = size;
    FormatMetadata mtdata = GetFormatMetadata<T>();
    new (&formats[index]) Format{.metadata = mtdata, .creator = &create_instance<T>, .is_valid = true};

    std::unordered_map<std::string, size_t>::const_iterator it;
    for (const std::string& alias : mtdata.aliases) {
      it = alias_to_index.find(alias);
      if (it != alias_to_index.end()) {
        continue; // alias already used TODO: Add warning message
      }
      alias_to_index.insert({alias, index});
    }
    ++size;
  }

  Format from_alias(const std::string& alias) const {
    std::unordered_map<std::string, size_t>::const_iterator it = alias_to_index.find(alias);
    if (it == alias_to_index.end())
      return Format{FormatMetadata{}, nullptr, false};
    return formats[it->second];
  }

  inline std::string infos() const {
    std::string infos = "Supported atom format:\n";
    for (size_t i = 0; i < size; ++i) {
      infos += strf("\t- %s\n", formats.at(i).metadata.name.c_str());
    }
    infos += "\n";
    return infos;
  }
};

} // namespace factory

template <typename T> const factory::IOFactory<T> io_factory() {
  static const factory::IOFactory<T> instance{};
  return instance;
}

/* ------------------------------------------------------------------------- */
/* PARSERS                                                                   */
/* ------------------------------------------------------------------------- */
namespace lammps_data {

using namespace enum_traits;
using namespace tokens;
using namespace tokens::numeric;

// using AtomDataParserFn = bool (*)(const TokenSet&, IOContext&, size_t);
// using KeywordParserFn = bool (*)(const TokenSet&, IOContext&);

enum class AtomStyle {
  Atomic,
  Charge,
  Molecular,
  Full,
  Count,
  None,
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

constexpr TokenNeedles<enum_size<AtomStyle>()> NEEDLES_ATOM_STYLE{
    "atomic",
    "charge",
    "molecule",
    "full",
};

// Needles for header keywords (HDRKEY)
constexpr TokenNeedles<1> NEEDLES_KEYWORD_ATOMS{"atoms"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_ATOM_TYPES{"atom", "types"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_XAXIS{"xlo", "xhi"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_YAXIS{"ylo", "yhi"};
constexpr TokenNeedles<2> NEEDLES_KEYWORD_ZAXIS{"zlo", "zhi"};
constexpr TokenNeedles<3> NEEDLES_KEYWORD_TILTS{"xy", "xz", "yz"};

constexpr TokenNeedles<enum_size<Section>()> NEEDLES_SECTIONS{
    "Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals",
};

template <EnumWithMemberCount Enum> struct ActiveKeywords {
  std::array<size_t, enum_size<Enum>()> indices{};
  size_t count{};

  constexpr ActiveKeywords() {
    for (size_t i = 0; i < enum_size<Keyword>(); ++i)
      indices[i] = i;
    count = enum_size<Keyword>();
  }

  inline void pop(size_t j) noexcept { indices[j] = indices[--count]; }
  inline bool empty() noexcept { return count == 0; }
};

template <size_t I, size_t J>
inline bool parse_keyword_box_axis(const TokenSet& tokens, IOContext& ctx, const TokenNeedles<2>& needles) {
  double boxlo, boxhi;
  if (!tokens.atleast(4))
    return false;
  if (!token_set_match_seq(tokens, needles))
    return false;
  if (!parse_single_token(tokens.at<0>(), boxlo))
    return false;
  if (!parse_single_token(tokens.at<1>(), boxhi))
    return false;
  matrix_traits::at<I, J>(ctx.data.cell) = boxhi - boxlo;
  return true;
}

template <Keyword K> inline bool parse_keyword(const TokenSet&, IOContext&) { return false; }

template <> inline bool parse_keyword<Keyword::Atoms>(const TokenSet& tokens, IOContext& ctx) {
  if (!tokens.atleast(2))
    return false;
  if (!token_match_any(tokens.at<1>(), NEEDLES_KEYWORD_ATOMS))
    return false;
  if (!parse_single_token(tokens.at<0>(), ctx.data.nat))
    return false;
  return true;
}

template <> inline bool parse_keyword<Keyword::AtomTypes>(const TokenSet& tokens, IOContext& ctx) {
  if (!tokens.atleast(3))
    return false;
  if (!token_set_match_seq(tokens, NEEDLES_KEYWORD_ATOM_TYPES))
    return false;
  if (!parse_single_token(tokens.at<0>(), ctx.data.ntypes))
    return false;
  return true;
}

template <> inline bool parse_keyword<Keyword::XAxis>(const TokenSet& tokens, IOContext& ctx) {
  return parse_keyword_box_axis<0, 0>(tokens, ctx, NEEDLES_KEYWORD_XAXIS);
}

template <> inline bool parse_keyword<Keyword::YAxis>(const TokenSet& tokens, IOContext& ctx) {
  return parse_keyword_box_axis<1, 1>(tokens, ctx, NEEDLES_KEYWORD_YAXIS);
}

template <> inline bool parse_keyword<Keyword::ZAxis>(const TokenSet& tokens, IOContext& ctx) {
  return parse_keyword_box_axis<2, 2>(tokens, ctx, NEEDLES_KEYWORD_ZAXIS);
}

template <> inline bool parse_keyword<Keyword::Tilts>(const TokenSet& tokens, IOContext& ctx) {
  double* xy = &ctx.data.cell.m12;
  double* xz = &ctx.data.cell.m13;
  double* yz = &ctx.data.cell.m23;

  if (!tokens.atleast(6))
    return false;
  if (!token_set_match_seq(tokens, NEEDLES_KEYWORD_TILTS))
    return false;
  if (!parse_single_token(tokens.at<0>(), *xy))
    return false;
  if (!parse_single_token(tokens.at<1>(), *xz))
    return false;
  if (!parse_single_token(tokens.at<2>(), *yz))
    return false;
  return true;
}

template <> inline bool parse_keyword<Keyword::Section>(const TokenSet& tokens, IOContext&) {
  if (tokens.atleast(1) && token_match_any(tokens.at<0>(), NEEDLES_SECTIONS))
    return true;
  return false;
}

template <size_t... Is> constexpr auto make_keyword_dispatch(std::index_sequence<Is...>) {
  return std::array{+[](const TokenSet& tokens, IOContext& ctx) {
    return parse_keyword<enum_from_index<Keyword>(Is)>(tokens, ctx);
  }...};
}

constexpr auto KeywordParsersDispatch = make_keyword_dispatch(std::make_index_sequence<enum_size<Keyword>()>{});

struct LAMMPSDataTraits : IOTraitsBase<LAMMPSDataTraits> {
  static constexpr const char* name = "LAMMPS Data";
  static constexpr std::array<const char*, 3> aliases{"data", "lmp-data", "lmp"};
  static constexpr const char* description = "LAMMPS Data format";
};

class LAMMPSDataParser : public TextParser {
public:
  using Traits = LAMMPSDataTraits;
  LAMMPSDataParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression), m_atom_style(AtomStyle::None), m_current_section(Section::Header) {}

  int64_t next() override {
    if (size_t cursor = m_file.tell(); cursor == 0) {
      m_file.get_line();
      return static_cast<int64_t>(cursor);
    } else {
      return -1;
    }
  }

  bool parse(IOContext& ctx) override {

    m_current_section = Section::Header;

    while (!m_file.eof()) {
      if (!skipline(current_line())) {
        switch (m_current_section) {
        case Section::Header:
          parse_section_header(ctx);
          if (ctx.has<IOContext::DOMAIN_ONLY>())
            return true;
          break;
        case Section::Masses:
          parse_section_masses();
          break;
        case Section::Atoms:
          parse_section_atoms(ctx);
          break;
        case Section::Velocities:
          parse_section_velocities(ctx);
          break;
        case Section::Ignored:
          jumpto_next_section();
          break;
        default:
          break;
        }
      }

      if (m_current_section != Section::Ignored) {
        current_line() = trim_spaces(m_file.get_line());
      }
    }
    return true;
  };

private:
  AtomStyle m_atom_style;
  Section m_current_section;

  // FIXME: The string_view point to internal file buffer, if buffer is flushed or moved it can lead
  // to invalid string_view
  inline std::string_view& section_line(void) { return m_line_buffer[1]; }
  inline const std::string_view& section_line(void) const { return m_line_buffer[1]; }

  bool parse_section_header(IOContext& ctx) {

    ActiveKeywords<Keyword> active{};

    while (!m_file.eof() && !active.empty()) {
      if (!skipline(current_line())) {
        tokenize(current_line(), m_tokens);

        for (size_t i = 0; i < active.count; ++i) {
          size_t idx = active.indices[i];
          if (KeywordParsersDispatch[idx](m_tokens, ctx)) {
            active.pop(i);

            if (idx == enum_to_index(Keyword::Section)) {
              m_current_section = Section::Ignored;
              return true;
            }

            break;
          }
        }
      }

      current_line() = trim_spaces(m_file.get_line());
    }

    m_current_section = Section::Ignored;
    return true;
  }

  template <bool RemapAtomID>
  inline bool parse_atom_line_atomic(const TokenSet& tokens, IOContext& ctx, size_t particle_count) {

    size_t index{};
    size_t id{}, type{};
    double x{}, y{}, z{};

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

    policy::Traits::assign_position(ctx.data, index, x, y, z);
    policy::Traits::assign_type(ctx.data, index, type - 1);

    return true;
  }

  template <bool RemapAtomID>
  inline bool parse_atom_line_charge(const TokenSet& tokens, IOContext& ctx, size_t particle_count) {

    size_t index{};
    size_t id{}, type{};
    double x{}, y{}, z{};

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

    policy::Traits::assign_position(ctx.data, index, x, y, z);
    policy::Traits::assign_type(ctx.data, index, type - 1);

    return true;
  }

  template <AtomStyle style, bool RemapAtomID> inline bool parse_section_atoms_impl(IOContext& ctx) {

    bool ok;
    size_t particle_count = 0;
    policy::Traits::allocate(ctx.data, ctx.data.nat);

    while (!m_file.eof()) {
      if (!skipline(current_line())) {
        tokenize(current_line(), m_tokens);

        if (parse_keyword<Keyword::Section>(m_tokens, ctx))
          break;

        bool atom_overflow = cmp::gt(++particle_count, ctx.data.nat);

        if constexpr (style == AtomStyle::Atomic) {
          ok = parse_atom_line_atomic<RemapAtomID>(m_tokens, ctx, particle_count);
        } else if constexpr (style == AtomStyle::Charge) {
          ok = parse_atom_line_charge<RemapAtomID>(m_tokens, ctx, particle_count);
        }

        // TODO: Change this
        if (atom_overflow || !ok) {
          if (atom_overflow) {
            PANIC("Too many atoms expected %ld got %ld", ctx.data.nat, particle_count);
          }
          PANIC("Error encounterd at line: %s", std::string(current_line()).c_str())
        }
      }

      current_line() = trim_spaces(m_file.get_line());
    }

    if (particle_count != ctx.data.nat) {
      PANIC("Number of particles is not consitent...");
    }

    m_current_section = Section::Ignored;
    return true;
  }

  bool parse_section_atoms(IOContext& ctx) {

    size_t index, style;
    tokenize(section_line(), m_tokens);
    if (token_set_match_any(m_tokens, NEEDLES_ATOM_STYLE, index, style)) {
      m_atom_style = enum_from_index<AtomStyle>(style);
    } else {
      PANIC("Invalid atom style defined.");
    }

    switch (m_atom_style) {
    case (AtomStyle::Atomic):
      if (ctx.has<IOContext::REMAP_ATOM>())
        return parse_section_atoms_impl<AtomStyle::Atomic, true>(ctx);
      else
        return parse_section_atoms_impl<AtomStyle::Atomic, false>(ctx);
      break;
    case (AtomStyle::Charge):
    default:
      PANIC("Not Implemented");
      break;
    }

    m_current_section = Section::Ignored;
    return false;
  }

  bool parse_section_masses(void) {
    m_current_section = Section::Ignored;
    return true;
  }

  template <bool RemapAtomID>
  inline bool parse_atom_velocities(const TokenSet& tokens, IOContext& ctx, size_t particle_count) {

    size_t id{}, index{};
    double vx{}, vy{}, vz{};

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

    policy::Traits::assign_velocity(ctx.data, index, vx, vy, vz);

    return true;
  }

  template <bool RemapAtomID> bool parse_section_velocities_impl(IOContext& ctx) {

    size_t particle_count = 0;

    while (!m_file.eof()) {
      if (!skipline(current_line())) {
        tokenize(current_line(), m_tokens);

        if (parse_keyword<Keyword::Section>(m_tokens, ctx))
          break;

        bool atom_overflow = cmp::gt(++particle_count, ctx.data.nat);

        if (atom_overflow || !parse_atom_velocities<RemapAtomID>(m_tokens, ctx, particle_count)) {
          if (atom_overflow) {
            PANIC("Too many atoms expected %ld got %ld", ctx.data.nat, particle_count);
          }
          PANIC("Error encounterd at line: %s", std::string(current_line()).c_str())
        }
      }

      current_line() = trim_spaces(m_file.get_line());
    }

    if (particle_count != ctx.data.nat) {
      PANIC("Number of particles is not consitent...");
    }

    m_current_section = Section::Ignored;
    return true;
  }

  bool parse_section_velocities(IOContext& ctx) {
    if (ctx.has<IOContext::REMAP_ATOM>())
      return parse_section_velocities_impl<true>(ctx);
    else
      return parse_section_velocities_impl<false>(ctx);
  }

  inline void jumpto_next_section(void) {
    current_line() = trim_spaces(current_line());
    while (!m_file.eof()) {
      if (!skipline(current_line())) {

        tokenize(current_line(), m_tokens);

        if (token_match_any(m_tokens[0], NEEDLES_SECTIONS)) {
          set_current_section(m_tokens[0]);
          section_line() = current_line(); // keep track of the section for later
          break;
        }
      }
      current_line() = trim_spaces(m_file.get_line());
    }
  }

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

namespace xyz {

using namespace enum_traits;
using namespace tokens;
using namespace tokens::numeric;
using namespace base_types;

enum class FieldType : uint8_t {
  R = 1 << 0,
  S = 1 << 1,
  I = 1 << 2,
  Count = 3,
  None = 0,
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

constexpr TokenNeedles<enum_size<FieldType>()> NEEDLES_FIELD_TYPE{"R", "S", "I"};
constexpr TokenNeedles<enum_size<Field>()> NEEDLES_FIELDS{
    "species", "pos", "velo", "id", "type",
};

struct Property {
  Field field = Field::None;
  size_t icol{};
  size_t ncol{};
};

struct XYZProperties : public properties_traits::GenericContainer<Property, Field, &Property::field> {
  using Enum = properties_traits::GenericContainer<Property, Field, &Property::field>::enum_t;

  template <Enum... args> constexpr void expect() {
    ((len += (has<args>()) ? get<args>().ncol : static_cast<size_t>(0)), ...);
  }

  IJK loc_positions{};
  IJK loc_velocities{};
};

enum class FieldFlag : uint16_t {
  ID = 1 << 0,
  Type = 1 << 1,
  Species = 1 << 2,
  Velocities = 1 << 3,
  RemapID = 1 << 4,
};

template <bool RemapAtomID, bool HasID, bool HasParticleSpecies, bool HasParticleType, bool HasVelocities>
struct LineParserTmpl {
  inline static bool call(const TokenSet& tokens, IOContext& ctx, const XYZProperties& ppt, size_t particle_count) {

    size_t id{}, index{};
    ssize_t type{};
    Vec3d r{}, v{};

    if (!tokens.atleast(ppt.len)) {
      linfo("ERR: line too shorts %ld %ld\n", tokens.size(), ppt.len);
      return false;
    }

    if constexpr (HasID) {
      if (!parse_single_token(tokens[ppt.get<Field::ID>().icol], id))
        return false;
    } else {
      id = particle_count;
    }

    const IJK& loc = ppt.loc_positions;
    if (!parse_single_token(tokens[loc.i], r.x))
      return false;
    if (!parse_single_token(tokens[loc.j], r.y))
      return false;
    if (!parse_single_token(tokens[loc.k], r.z))
      return false;

    if constexpr (HasParticleSpecies) {
      type = 0;
    } else if constexpr (HasParticleType) {
      if (!parse_single_token(tokens[ppt.get<Field::Type>().icol], type))
        return false;
    } else {
      type = 0;
    }

    if constexpr (HasVelocities) {
      const IJK& vloc = ppt.loc_velocities;
      if (!parse_single_token(tokens[vloc.i], v.x))
        return false;
      if (!parse_single_token(tokens[vloc.j], v.y))
        return false;
      if (!parse_single_token(tokens[vloc.k], v.z))
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

    policy::Traits::assign_position(ctx.data, index, r.x, r.y, r.z);
    policy::Traits::assign_type(ctx.data, index, type);

    if constexpr (HasVelocities) {
      policy::Traits::assign_velocity(ctx.data, index, v.x, v.y, v.z);
    }

    return true;
  }
};

template <uint16_t Mask> constexpr auto make_parser_dispatch() {

  constexpr bool has_id = Mask & to_underlying(FieldFlag::ID);
  constexpr bool has_species = Mask & to_underlying(FieldFlag::Species);
  constexpr bool has_type = Mask & to_underlying(FieldFlag::Type);
  constexpr bool remap = Mask & to_underlying(FieldFlag::RemapID);
  constexpr bool has_velocities = Mask & to_underlying(FieldFlag::Velocities);

  using LineParser = LineParserTmpl<remap, has_id, has_species, has_type, has_velocities>;

  return +[](const TokenSet& tokens, IOContext& ctx, const XYZProperties& ppt, size_t pc) -> bool {
    return LineParser::call(tokens, ctx, ppt, pc);
  };
}

template <size_t... Is> constexpr auto build_dispatch_table_impl(std::index_sequence<Is...>) {
  return std::array{make_parser_dispatch<Is>()...};
}

template <size_t N> constexpr auto build_dispatch_table() {
  return build_dispatch_table_impl(std::make_index_sequence<(1 << N)>{});
}

constexpr auto parser_dispatch = build_dispatch_table<5>();

template <size_t N, size_t M>
bool parse_extxyz_comment_line(StringView line, TokenSetTmpl<N>& tok1, TokenSetTmpl<M>& tok2) {

  constexpr StringView key_lattice = "Lattice=\"";
  constexpr size_t key_lattice_len = key_lattice.size(); // 9
  constexpr StringView key_properties = "Properties=";
  constexpr size_t key_properties_len = key_properties.size(); // 11

  const char* begin = line.data();
  const char* end = line.data() + line.size();

  if (!startswith(line, key_lattice))
    return false;

  // 1. Parse Lattice=...
  const char* first = begin + key_lattice_len;
  const char* last = first;

  while (last < end && !DoubleQuoteDelimiter{}(last))
    ++last;

  if (last >= end)
    return false;

  tokenize(StringView(first, last), tok1);

  // 2. Parse Properties=....
  first = last;
  while (first + key_properties_len <= end) {

    if (std::memcmp(first, key_properties.data(), key_properties_len) == 0) {
      first += key_properties_len;
      last = first;

      while (last < end && !SpaceDelimiter{}(last)) {
        ++last;
      }

      tokenize_impl(StringView(first, last), tok2, ColumnDelimiter{});
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

inline bool parse_extxyz_properties(const TokenSet& tokens, XYZProperties& properties, FieldFlag& flags) {

  if (tokens.size() % 3 != 0)
    return false;

  size_t icol = 0;
  size_t ncol = 0;
  size_t ifield = 0;
  size_t itype = 0;

  properties.clear();

  for (size_t i = 0; i < tokens.size(); i += 3) {
    size_t j = i + 1;
    size_t k = i + 2;

    parse_single_token(tokens[k], ncol);

    if ((token_match_any(tokens[i], NEEDLES_FIELDS, ifield) && token_match_any(tokens[j], NEEDLES_FIELD_TYPE, itype))) {
      properties.set({
          .field = enum_from_index<Field>(ifield),
          .icol = icol,
          .ncol = ncol,
      });
    }

    icol += ncol;
  }

  if (!properties.has<Field::Positions>())
    PANIC("Positions are not defined in atom properties");

  ssize_t ipos = static_cast<ssize_t>(properties.get<Field::Positions>().icol);
  properties.loc_positions = {ipos, ipos + 1, ipos + 2};

  if (properties.has<Field::ID>()) {
    flags |= FieldFlag::ID;
    properties.expect<Field::ID>();
  }

  if (!(properties.has<Field::Type>() || properties.has<Field::Species>())) {
    PANIC("Either Species or TYPE should be define in xyz properties")
  } else if (properties.has<Field::Species>()) {
    flags |= FieldFlag::Species;
    properties.expect<Field::Species>();
  } else if (properties.has<Field::Type>()) {
    flags |= FieldFlag::Type;
    properties.expect<Field::Type>();
  }

  if (properties.has<Field::Velocities>()) {
    flags |= FieldFlag::Velocities;
    ssize_t ivel = static_cast<ssize_t>(properties.get<Field::Velocities>().icol);
    properties.loc_velocities = {ivel, ivel + 1, ivel + 2};
    properties.expect<Field::Velocities>();
  }

  return true;
}

struct ExtendedXYZTraits : IOTraitsBase<ExtendedXYZTraits> {
  static constexpr const char* name = "Extended XYZ";
  static constexpr std::array<const char*, 3> aliases{"xyz", "exyz", "ext-xyz"};
  static constexpr const char* description = "Extended XYZ";
};

class ExtendedXYZParser : public TextParser {
public:
  using Traits = ExtendedXYZTraits;
  ExtendedXYZParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override {
    if (size_t cursor = m_file.tell(); cursor == 0) {
      m_file.get_line();
      return static_cast<int64_t>(cursor);
    } else {
      return -1;
    }
  }

  bool parse(IOContext& ctx) override {

    current_line() = {};
    while (!m_file.eof() && skipline(current_line()))
      current_line() = trim_spaces(m_file.get_line());

    // check if the first non empty line is a valid number
    tokenize(current_line(), m_tokens);
    if (!parse_tokens(m_tokens.at<0>(), ctx.data.nat)) {
      PANIC("Invalid atom number: '%s'", std::string(current_line()).c_str());
    }

    // next line should be the comment line with lattice and properties definitions
    current_line() = trim_spaces(m_file.get_line());
    if (!parse_extxyz_comment_line(current_line(), m_tokens, m_ptokens))
      PANIC("Failed to read the comment line: '%s'", std::string(current_line()).c_str());

    if (!parse_extxyz_lattice(m_tokens, ctx))
      PANIC("Failed to parse lattice.");

    if (ctx.has<IOContext::DOMAIN_ONLY>())
      return true;

    FieldFlag mask = static_cast<FieldFlag>(0);
    if (ctx.has<IOContext::REMAP_ATOM>())
      mask |= FieldFlag::RemapID;

    if (!parse_extxyz_properties(m_ptokens, m_properties, mask))
      PANIC("Filed to parse properties.")

    size_t particle_count = 0;
    auto& parse_line = parser_dispatch[to_underlying(mask)];
    policy::Traits::allocate(ctx.data, ctx.data.nat);

    current_line() = trim_spaces(m_file.get_line());
    while (!m_file.eof() && particle_count < ctx.data.nat) {
      if (!skipline(current_line())) {
        tokenize(current_line(), m_tokens);

        bool atom_overflow = cmp::gt(++particle_count, ctx.data.nat);
        if (atom_overflow || !parse_line(m_tokens, ctx, m_properties, particle_count)) {
          if (atom_overflow)
            PANIC("Too many atoms expected %ld got %ld", ctx.data.nat, particle_count);
          PANIC("Error encounterd at line: %s", std::string(current_line()).c_str())
        }
      }
      current_line() = trim_spaces(m_file.get_line());
    }

    if (particle_count != ctx.data.nat)
      PANIC("Number of particles is not consitent...");

    return true;
  }

private:
  TokenSet m_ptokens; // token for property parsing
  XYZProperties m_properties;
};

} // namespace xyz

namespace lammps_dump {

using namespace base_types;
using namespace tokens;
using namespace tokens::numeric;
using namespace enum_traits;

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
    "id", "mol", "type", "element", "x",   "y",   "z",  "xu", "yu", "zu",
    "xs", "ys",  "zs",   "xsu",     "ysu", "zsu", "vx", "vy", "vz",
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

  ctx.data.cell.m21 = 0.0;
  ctx.data.cell.m31 = 0.0;
  ctx.data.cell.m32 = 0.0;

  ctx.data.origin.x = xlo;
  ctx.data.origin.y = ylo;
  ctx.data.origin.z = zlo;
};

struct Property {
  Field field = Field::None;
  size_t index{};
};

struct Properties : public properties_traits::GenericContainer<Property, Field, &Property::field> {
  using Enum = properties_traits::GenericContainer<Property, Field, &Property::field>::enum_t;

  Coordinates coord{};
  IJK rloc{}, vloc{};

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

enum class FieldFlag : uint16_t {
  Triclinic = 1 << 0,
  Unwrapped = 1 << 1,
  Scaled = 1 << 2,
  ID = 1 << 3,
  MolID = 1 << 4,
  Type = 1 << 5,
  Velocities = 1 << 6,
  RemapAtom = 1 << 7,
};

inline bool parse_atom_properties(const TokenSet& tokens, Properties& properties, FieldFlag& mask) {

  constexpr size_t headlen = 2; // ITEM: ATOMS
  size_t ifield;

  properties.clear();

  for (size_t i = headlen; i < tokens.size(); ++i) {
    if (token_match_any(tokens[i], NEEDLES_FIELDS, ifield)) {
      properties.set({.field = enum_from_index<Field>(ifield), .index = i - headlen});
    }
  }

  if (properties.has<Field::ID>())
    mask |= FieldFlag::ID;

  if (properties.has<Field::MolID>())
    mask |= FieldFlag::MolID;

  if (properties.has<Field::Type>())
    mask |= FieldFlag::Type;

  if (!properties.has_positions()) {
    PANIC("No positions fields are present...");
  }

  // deduced coodinate system base on present field
  for (const Coordinates& c : COORDINATES) {
    if (properties.has(c.x, c.y, c.z)) {
      properties.coord = c;
      properties.rloc = properties.ijk(c.x, c.y, c.z);
      properties.expect(c.x, c.y, c.z);
      break;
    }
  }

  if (properties.coord.scaled)
    mask |= FieldFlag::Scaled;
  if (properties.coord.unwrapped)
    mask |= FieldFlag::Unwrapped;

  // check if there are velocity fields
  if (properties.has<Field::VX, Field::VY, Field::VZ>()) {
    mask |= FieldFlag::Velocities;
    properties.vloc = properties.ijk<Field::VX, Field::VY, Field::VZ>();
    properties.expect(Field::VX, Field::VY, Field::VZ);
  }

  return true;
}

template <bool IsTriclinic, matrix_traits::Vec3dLike Vec> struct BBoxParserTmpl {
  inline bool operator()(const TokenSet& tokens, const Vec& v) {
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
      v.z = 0.0;
      return true;
    }
  }
};

template <bool IsTriclinic, bool Scaled, bool Unwrapped, bool HasID, bool HasMolID, bool HasType, bool HasVelocities,
          bool RemapAtomID>
struct ParserTmpl {
  inline static bool call(const TokenSet& tokens, IOContext& ctx, Properties& ppt, size_t particle_count) {
    size_t id{}, index{};
    ssize_t type{};
    Vec3d r{}, v{};

    if (!tokens.atleast(ppt.len)) {
      linfo("ERR: line too shorts %ld %ld\n", tokens.size(), ppt.len);
      return false;
    }

    if constexpr (HasID) {
      if (!parse_single_token(tokens[ppt.get<Field::ID>().index], id))
        return false;
    } else {
      id = particle_count;
    }

    if constexpr (HasType) {
      if (!parse_single_token(tokens[ppt.get<Field::Type>().index], type))
        return false;
    } else {
      type = 0;
    }

    if (!parse_single_token(tokens[ppt.rloc.i], r.x))
      return false;
    if (!parse_single_token(tokens[ppt.rloc.j], r.y))
      return false;
    if (!parse_single_token(tokens[ppt.rloc.k], r.z))
      return false;

    if constexpr (Scaled) {
      if constexpr (IsTriclinic) {
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

    if constexpr (HasVelocities) {
      if (!parse_single_token(tokens[ppt.vloc.i], v.x))
        return false;
      if (!parse_single_token(tokens[ppt.vloc.j], v.y))
        return false;
      if (!parse_single_token(tokens[ppt.vloc.k], v.z))
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

    policy::Traits::assign_position(ctx.data, index, r.x, r.y, r.z);
    policy::Traits::assign_type(ctx.data, index, type);

    if constexpr (HasVelocities) {
      policy::Traits::assign_velocity(ctx.data, index, v.x, v.y, v.z);
    }

    return true;
  }
};

template <uint16_t Mask> constexpr auto make_parser_dispatch() {
  constexpr bool is_triclinic = Mask & enum_traits::to_underlying(FieldFlag::Triclinic);
  constexpr bool is_scaled = Mask & enum_traits::to_underlying(FieldFlag::Scaled);
  constexpr bool is_unwrapped = Mask & enum_traits::to_underlying(FieldFlag::Unwrapped);
  constexpr bool has_id = Mask & enum_traits::to_underlying(FieldFlag::ID);
  constexpr bool has_molid = Mask & enum_traits::to_underlying(FieldFlag::MolID);
  constexpr bool has_type = Mask & enum_traits::to_underlying(FieldFlag::Type);
  constexpr bool has_velocities = Mask & enum_traits::to_underlying(FieldFlag::Velocities);
  constexpr bool remap_atom = Mask & enum_traits::to_underlying(FieldFlag::RemapAtom);

  using Parser =
      ParserTmpl<is_triclinic, is_scaled, is_unwrapped, has_id, has_molid, has_type, has_velocities, remap_atom>;
  return +[](const TokenSet& tokens, IOContext& ctx, Properties& ppt, size_t pc) -> bool {
    return Parser::call(tokens, ctx, ppt, pc);
  };
}

template <size_t... Is> constexpr auto build_dispatch_table_impl(std::index_sequence<Is...>) {
  return std::array{make_parser_dispatch<Is>()...};
}

template <size_t N> constexpr auto build_dispatch_table() {
  return build_dispatch_table_impl(std::make_index_sequence<(1 << N)>{});
}

constexpr auto parser_dispatch = build_dispatch_table<8>();

struct LAMMPSDumpTraits : IOTraitsBase<LAMMPSDumpTraits> {
  static constexpr const char* name = "LAMMPS Dump";
  static constexpr std::array<const char*, 3> aliases{"dump", "lmp-dump", "lammps-dump"};
  static constexpr const char* description = "LAMMPS Dump format";
};

class LAMMPSDumpParser : public TextParser {
public:
  using Traits = LAMMPSDumpTraits;
  LAMMPSDumpParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override {
    if (size_t cursor = m_file.tell(); cursor == 0) {
      m_file.get_line();
      return static_cast<int64_t>(cursor);
    } else {
      return -1;
    }
  }

  inline void skip_to_label(const StringView str) {
    while (!m_file.eof()) {
      if (startswith(current_line(), str))
        break;
      current_line() = trim_spaces(m_file.get_line());
    }
  }

  bool parse(IOContext& ctx) override {
    constexpr StringView label_natom = "ITEM: NUMBER OF ATOMS";
    constexpr StringView label_box_bounds = "ITEM: BOX BOUNDS";
    constexpr StringView label_atoms = "ITEM: ATOMS";

    current_line() = {};

    Properties properties;
    FieldFlag mask = static_cast<FieldFlag>(0);

    if (ctx.has<IOContext::REMAP_ATOM>())
      mask |= FieldFlag::RemapAtom;

    skip_to_label(label_natom);

    // line just after 'ITER: NUMBER OF ATOMS' should be a valid int
    current_line() = trim_spaces(m_file.get_line());
    tokenize(current_line(), m_tokens, 1);
    if (!parse_single_token(m_tokens.at<0>(), ctx.data.nat)) {
      PANIC("Invalid number of atom: %s", std::string(current_line()).c_str());
    }

    // 'ITEM: BOX BOUNDS' is expected after NUMBER OF ATOMS
    skip_to_label(label_box_bounds);
    if (!startswith(current_line(), label_box_bounds))
      PANIC("INVALID BOX BOUNDS: '%s'", std::string(current_line()).c_str());

    tokenize(current_line(), m_tokens);
    ctx.set<IOContext::TRICLINIC>(m_tokens.atleast(7));

    Mat3d bbox{};

    if (ctx.has<IOContext::TRICLINIC>()) {
      mask |= FieldFlag::Triclinic;
      auto parser = BBoxParserTmpl<true, matrix_traits::Vec3dView>{};
      for (size_t i = 0; i < 3; ++i) {
        tokenize(trim_spaces(m_file.get_line()), m_tokens);
        if (!parser(m_tokens, matrix_traits::row(i, bbox)))
          PANIC("Faild to parse 'ATOM: BOX BOUNDS' at line %*s", current_line().size(), current_line().data());
      }
    } else {
      auto parser = BBoxParserTmpl<false, matrix_traits::Vec3dView>{};
      for (size_t i = 0; i < 3; ++i) {
        tokenize(trim_spaces(m_file.get_line()), m_tokens);
        if (!parser(m_tokens, matrix_traits::row(i, bbox)))
          PANIC("Faild to parse 'ATOM: BOX BOUNDS' at line %*s", current_line().size(), current_line().data());
      }
    }

    convert_bbox_to_lattice_vectors(ctx, bbox);
    // early return if we only care about the lattice vector
    if (ctx.has<IOContext::DOMAIN_ONLY>())
      return true;

    // 'ITEM: ATOMS' is expected after BOX BOUNDS
    skip_to_label(label_atoms);
    if (!startswith(current_line(), label_atoms))
      PANIC("INVALID ATOMS: '%s'", std::string(current_line()).c_str());

    tokenize(current_line(), m_tokens);
    parse_atom_properties(m_tokens, properties, mask);

    size_t particle_count = 0;
    auto& parse_line = parser_dispatch[to_underlying(mask)];
    policy::Traits::allocate(ctx.data, ctx.data.nat);

    current_line() = trim_spaces(m_file.get_line());
    while (!m_file.eof() && particle_count < ctx.data.nat) {
      if (!skipline(current_line())) {
        tokenize(current_line(), m_tokens);

        bool atom_overflow = cmp::gt(++particle_count, ctx.data.nat);
        if (atom_overflow || !parse_line(m_tokens, ctx, properties, particle_count)) {
          if (atom_overflow)
            PANIC("Too many atoms expected %ld got %ld", ctx.data.nat, particle_count);
          PANIC("Error encounterd at line: %s", std::string(current_line()).c_str())
        }
      }
      current_line() = trim_spaces(m_file.get_line());
    }

    return true;
  }
};
} // namespace lammps_dump

namespace vasp {

using namespace base_types;
using namespace tokens;
using namespace tokens::numeric;
using namespace enum_traits;

constexpr TokenNeedles<1> CARTESIAN_COORDINATES{"Cartesian"};

enum class FieldFlag : uint16_t {
  Cartesian = 1 << 0,
};

template <bool IsCartesian> struct LineParserTmpl {
  inline static bool call(const TokenSet& tokens, IOContext& ctx, size_t particle_count, const Vec3d& scale) {

    Vec3d r{};
    size_t index = particle_count - 1;

    if (!tokens.atleast(3)) {
      linfo("ERR: line too shorts %ld %ld\n", tokens.size(), 3);
      return false;
    }

    if (!parse_single_token(tokens.at<0>(), r.x))
      return false;
    if (!parse_single_token(tokens.at<1>(), r.y))
      return false;
    if (!parse_single_token(tokens.at<2>(), r.z))
      return false;

    if constexpr (IsCartesian) {
      r.x = r.x * scale.x;
      r.y = r.y * scale.y;
      r.z = r.z * scale.z;
    } else {
      r.x = r.x * ctx.data.cell.m11 + r.y * ctx.data.cell.m12 + r.z * ctx.data.cell.m13;
      r.x = r.x * ctx.data.cell.m21 + r.y * ctx.data.cell.m22 + r.z * ctx.data.cell.m23;
      r.x = r.x * ctx.data.cell.m31 + r.y * ctx.data.cell.m32 + r.z * ctx.data.cell.m33;
    }

    size_t type = ctx.species.type_from_index(index);

    policy::Traits::assign_position(ctx.data, index, r.x, r.y, r.z);
    policy::Traits::assign_type(ctx.data, index, type);

    return true;
  }
};

template <uint16_t Mask> constexpr auto make_parser_dispatch() {
  constexpr bool cartesian = Mask & to_underlying(FieldFlag::Cartesian);
  using LineParser = LineParserTmpl<cartesian>;
  return +[](const TokenSet& tokens, IOContext& ctx, size_t pc, const Vec3d& scaling) -> bool {
    return LineParser::call(tokens, ctx, pc, scaling);
  };
}

template <size_t... Is> constexpr auto build_dispatch_table_impl(std::index_sequence<Is...>) {
  return std::array{make_parser_dispatch<Is>()...};
}

template <size_t N> constexpr auto build_dispatch_table() {
  return build_dispatch_table_impl(std::make_index_sequence<(1 << N)>{});
}

constexpr auto parser_dispatch = build_dispatch_table<1>();

template <bool Milady> struct VaspTraits : IOTraitsBase<VaspTraits<Milady>> {};

template <> struct VaspTraits<false> : IOTraitsBase<VaspTraits<false>> {
  static constexpr const char* name = "VASP POSCAR";
  static constexpr std::array<const char*, 3> aliases{"vasp", "poscar", "pos"};
  static constexpr const char* description = "???";
};

template <> struct VaspTraits<true> : IOTraitsBase<VaspTraits<true>> {
  static constexpr const char* name = "VASP POSCAR MILADY";
  static constexpr std::array<const char*, 2> aliases{"mld", "milady"};
  static constexpr const char* description = "???";
};

template <bool Milady> class VaspPoscarParserTmpl : public TextParser {
public:
  using Traits = VaspTraits<Milady>;
  VaspPoscarParserTmpl(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override {
    if (size_t cursor = m_file.tell(); cursor == 0) {
      m_file.get_line();
      return static_cast<int64_t>(cursor);
    } else {
      return -1;
    }
  }

  bool parse(IOContext& ctx) override {

    // alway skip the first line, this is a comment
    m_file.get_line();

    // go to the first non empty line
    while (!m_file.eof() && skipline(current_line())) {
      current_line() = trim_spaces(m_file.get_line());
    }

    // parse scaling factor
    Vec3d scaling{};
    tokenize(current_line(), m_tokens);
    if (!(m_tokens.size() == 1 || m_tokens.size() == 3)) {
      PANIC("Ill formating sclaing factor");
    }

    if (m_tokens.size() == 1) {
      if (!parse_single_token(m_tokens.at<0>(), scaling.x))
        return false;
      scaling.y = scaling.x;
      scaling.z = scaling.x;
    } else {
      if (!parse_single_token(m_tokens.at<0>(), scaling.x))
        return false;
      if (!parse_single_token(m_tokens.at<1>(), scaling.y))
        return false;
      if (!parse_single_token(m_tokens.at<2>(), scaling.z))
        return false;
    }

    // parse lattice vectors
    // ax ay az
    // bx by bz
    // cx cy cz
    // /!\ a, b and c are stored column wise. H = (a | b | c)
    for (size_t i = 0; i < 3; ++i) {
      tokenize(trim_spaces(m_file.get_line()), m_tokens);
      matrix_traits::Vec3dView latvec = matrix_traits::column(i, ctx.data.cell);
      if (!m_tokens.atleast(3))
        return false;
      if (!parse_single_token(m_tokens.at<0>(), latvec.x))
        return false;
      if (!parse_single_token(m_tokens.at<1>(), latvec.y))
        return false;
      if (!parse_single_token(m_tokens.at<2>(), latvec.z))
        return false;
      matrix_traits::scale(latvec, matrix_traits::at(scaling, i));
    }

    // parse species: Al Cu ...
    tokenize(trim_spaces(m_file.get_line()), m_tokens);
    if (!m_tokens.atleast(1))
      return false;
    for (size_t i = 0; i < m_tokens.size(); ++i)
      ctx.species.add_species(m_tokens[i], i);

    // count per species
    tokenize(trim_spaces(m_file.get_line()), m_tokens);
    size_t count{};
    if (!m_tokens.atleast(ctx.species.size()))
      return false;
    for (size_t i = 0; i < ctx.species.size(); ++i) {
      if (!parse_single_token(m_tokens[i], count))
        return false;
      ctx.data.nat += count;
      ctx.data.ntypes += 1;
      ctx.species.set_type_count(i, count);
    }

    ctx.species.compute_count_offsets();

    // Parse coordinates system, Cartesian/Direct
    tokenize(trim_spaces(m_file.get_line()), m_tokens);
    bool xcart = token_match_any_no_case(m_tokens[0], CARTESIAN_COORDINATES);

    FieldFlag mask = static_cast<FieldFlag>(0);
    if (xcart)
      mask |= FieldFlag::Cartesian;

    auto& parse_line = parser_dispatch[to_underlying(mask)];
    size_t particle_count = 0;
    policy::Traits::allocate(ctx.data, ctx.data.nat);

    current_line() = trim_spaces(m_file.get_line());
    while (!m_file.eof() && particle_count < ctx.data.nat) {
      if (!skipline(current_line())) {
        tokenize(current_line(), m_tokens);

        bool atom_overflow = cmp::gt(++particle_count, ctx.data.nat);
        if (atom_overflow || !parse_line(m_tokens, ctx, particle_count, scaling)) {
          if (atom_overflow)
            PANIC("Too many atoms expected %ld got %ld", ctx.data.nat, particle_count);
          PANIC("Error encounterd at line: %s", std::string(current_line()).c_str())
        }
      }
      current_line() = trim_spaces(m_file.get_line());
    }

    if (particle_count != ctx.data.nat)
      PANIC("Number of particles is not consitent...");

    return true;
  }
};

template <bool Milady> class VaspPoscarWriterTmpl : public TextWriter {
public:
  using Traits = VaspTraits<Milady>;
  VaspPoscarWriterTmpl(std::string filepath, FileMode mode, FileCompression compression)
      : TextWriter(filepath, mode, compression) {}

  bool write(IOContext& ctx) override {
    m_buf.format = write_utils::FloatFormat{.width = 15, .fill = ' '};

    // comment line
    if constexpr (Milady) {
      m_buf.insert("111 ");
      m_buf.insert(ctx.species.size(), {.width = 0});
      for (size_t i = 0; i < ctx.species.size(); ++i) {
        m_buf.space();
        m_buf.insert(ctx.species.get_symbol(i));
        m_buf.space();
        m_buf.insert(0, {.width = 0});
      }
      m_buf.insert(" 0.00 0.00 0.00\n");
    } else {
      m_buf.insert("Auto-generated POSCAR file\n");
    }

    // sclaing factor
    m_buf.insert(1.0).newline();

    // lattice vector
    m_buf.insert(matrix_traits::column(0, ctx.data.cell)).newline();
    m_buf.insert(matrix_traits::column(1, ctx.data.cell)).newline();
    m_buf.insert(matrix_traits::column(2, ctx.data.cell)).newline();

    // Species
    for (size_t i = 0; i < ctx.species.size(); ++i) {
      m_buf.space();
      m_buf.insert(ctx.species.get_symbol(i));
    }
    m_buf.newline();

    // species count
    for (size_t i = 0; i < ctx.species.size(); ++i) {
      m_buf.space();
      m_buf.insert(ctx.species.get_count(i), {});
    }
    m_buf.newline();

    // FIXME: Ptr access with exastamp policy
    // write the positions
    // m_buf.insert("Cartesian").newline();
    // double* __restrict ptr = ctx.data.positions.data();
    // for (size_t i = 0; i < ctx.data.nat; ++i) {
    //   m_buf.insert(ptr[i * 3 + 0]).space();
    //   m_buf.insert(ptr[i * 3 + 1]).space();
    //   m_buf.insert(ptr[i * 3 + 1]).newline();
    // }

    if constexpr (Milady) {
      // Forces
      m_buf.newline();
      Vec3d f{};
      for (size_t i = 0; i < ctx.data.nat; ++i) {
        m_buf.insert(f).newline();
      }

      // Stress and spin flag
      m_buf.newline();
      m_buf.insert("0.0 0.0 0.0 0.0 0.0 0.0\n\n0\n");
    }

    m_buf.flush();
    return true;
  }
};

using VaspPoscarParser = VaspPoscarParserTmpl<false>;
using VaspMiladyParser = VaspPoscarParserTmpl<true>;

using VaspPoscarWriter = VaspPoscarWriterTmpl<false>;
using VaspMiladyWriter = VaspPoscarWriterTmpl<true>;

} // namespace vasp

namespace atomeye {

using namespace base_types;
using namespace tokens;
using namespace tokens::numeric;
using namespace enum_traits;

enum class Keyword {
  NumberOfParticles,
  LatticeH11,
  LatticeH12,
  LatticeH13,
  LatticeH21,
  LatticeH22,
  LatticeH23,
  LatticeH31,
  LatticeH32,
  LatticeH33,
  VelocityFlag,
  EntryCount,
  Count,
};

constexpr TokenNeedles<1> NEEDLES_NUM_OF_PARTICLES{"Number of particles"};
constexpr TokenNeedles<1> NEEDLES_LATTICE_H11{"H0(1,1)"};
constexpr TokenNeedles<1> NEEDLES_LATTICE_H12{"H0(1,2)"};
constexpr TokenNeedles<1> NEEDLES_LATTICE_H13{"H0(1,3)"};
constexpr TokenNeedles<1> NEEDLES_LATTICE_H21{"H0(2,1)"};
constexpr TokenNeedles<1> NEEDLES_LATTICE_H22{"H0(2,2)"};
constexpr TokenNeedles<1> NEEDLES_LATTICE_H23{"H0(2,3)"};
constexpr TokenNeedles<1> NEEDLES_LATTICE_H31{"H0(3,1)"};
constexpr TokenNeedles<1> NEEDLES_LATTICE_H32{"H0(3,2)"};
constexpr TokenNeedles<1> NEEDLES_LATTICE_H33{"H0(3,3)"};
constexpr TokenNeedles<1> NEEDLES_NO_VELOCITY{".NO_VELOCITY."};
constexpr TokenNeedles<1> NEEDLES_ENTRY_COUNT{"entry_count"};

template <EnumWithMemberCount Enum> struct ActiveKeywords {
  std::array<size_t, enum_size<Enum>()> indices{};
  size_t count{};

  constexpr ActiveKeywords() {
    for (size_t i = 0; i < enum_size<Keyword>(); ++i)
      indices[i] = i;
    count = enum_size<Keyword>();
  }

  inline void pop(size_t j) noexcept { indices[j] = indices[--count]; }
  inline bool empty() noexcept { return count == 0; }
};

template <size_t I, size_t J>
inline bool parse_keyword_h0(const TokenSet& tokens, IOContext& ctx, const TokenNeedles<1>& needles) {
  if (!tokens.atleast(2))
    return false;
  if (!token_match_any(tokens.at<0>(), needles))
    return false;
  if (!parse_single_token(cut_at_first(tokens.at<1>(), IsChar<' '>{}), matrix_traits::at<I, J>(ctx.data.cell)))
    return false;
  return true;
}

template <Keyword K> inline bool parse_keyword(const TokenSet&, IOContext&) { return false; }

template <> inline bool parse_keyword<Keyword::NumberOfParticles>(const TokenSet& tokens, IOContext& ctx) {
  if (!tokens.atleast(2))
    return false;
  if (!token_match_any(tokens.at<0>(), NEEDLES_NUM_OF_PARTICLES))
    return false;
  if (!parse_single_token(tokens.at<1>(), ctx.data.nat))
    return false;
  return true;
}

#define DEFINE_LATTICE_KEYWORD(K, I, J, NEEDLE)                                                                        \
  template <> inline bool parse_keyword<Keyword::K>(const TokenSet& tokens, IOContext& ctx) {                          \
    return parse_keyword_h0<I, J>(tokens, ctx, NEEDLE);                                                                \
  }

DEFINE_LATTICE_KEYWORD(LatticeH11, 0, 0, NEEDLES_LATTICE_H11);
DEFINE_LATTICE_KEYWORD(LatticeH12, 1, 0, NEEDLES_LATTICE_H12);
DEFINE_LATTICE_KEYWORD(LatticeH13, 2, 0, NEEDLES_LATTICE_H13);
DEFINE_LATTICE_KEYWORD(LatticeH21, 0, 1, NEEDLES_LATTICE_H21);
DEFINE_LATTICE_KEYWORD(LatticeH22, 1, 1, NEEDLES_LATTICE_H22);
DEFINE_LATTICE_KEYWORD(LatticeH23, 2, 1, NEEDLES_LATTICE_H23);
DEFINE_LATTICE_KEYWORD(LatticeH31, 0, 2, NEEDLES_LATTICE_H31);
DEFINE_LATTICE_KEYWORD(LatticeH32, 1, 2, NEEDLES_LATTICE_H32);
DEFINE_LATTICE_KEYWORD(LatticeH33, 2, 2, NEEDLES_LATTICE_H33);

#undef DEFINE_LATTICE_KEYWORD

template <> inline bool parse_keyword<Keyword::VelocityFlag>(const TokenSet& tokens, IOContext& ctx) {
  if (token_match_any(tokens.at<0>(), NEEDLES_NO_VELOCITY)) {
    ctx.set<IOContext::NO_VELOCITY>(true);
    return true;
  }
  return false;
}

template <> inline bool parse_keyword<Keyword::EntryCount>(const TokenSet& tokens, IOContext& ctx) {
  if (!tokens.atleast(2))
    return false;
  if (!token_match_any(tokens.at<0>(), NEEDLES_ENTRY_COUNT))
    return false;
  if (!parse_single_token(tokens.at<1>(), ctx.data.ntypes))
    return false;
  return true;
}

template <size_t... Is> constexpr auto make_keyword_dispatch(std::index_sequence<Is...>) {
  return std::array{+[](const TokenSet& tokens, IOContext& ctx) {
    return parse_keyword<enum_from_index<Keyword>(Is)>(tokens, ctx);
  }...};
}

constexpr auto KeywordParsersDispatch = make_keyword_dispatch(std::make_index_sequence<enum_size<Keyword>()>{});

struct AtomEyeTraits : IOTraitsBase<AtomEyeTraits> {
  static constexpr const char* name = "VASP POSCAR MILADY";
  static constexpr std::array<const char*, 1> aliases{"cfg"};
  static constexpr const char* description = "???";
};

class AtomEyeParser : public TextParser {
public:
  using Traits = AtomEyeTraits;
  AtomEyeParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override {
    if (size_t cursor = m_file.tell(); cursor == 0) {
      m_file.get_line();
      return static_cast<int64_t>(cursor);
    } else {
      return -1;
    }
  }

  bool parse_header(IOContext& ctx, size_t& entry_count) {
    ActiveKeywords<Keyword> active{};
    while (!m_file.eof() && !active.empty()) {
      if (!skipline(current_line())) {
        tokenize_impl(current_line(), m_tokens, EqualSignDelimiter{}, 3);
        for (size_t i = 0; i < active.count; ++i) {
          size_t idx = active.indices[i];
          if (KeywordParsersDispatch[idx](m_tokens, ctx)) {
            active.pop(i);

            if (idx == enum_to_index(Keyword::EntryCount)) {
              entry_count = ctx.data.ntypes;
              ctx.data.ntypes = 0;
              if (ctx.data.nat == 0) {
                linfo("Invalid number of particles");
                return false;
              }
              if (ctx.has<IOContext::NO_VELOCITY>() && entry_count < 3) {
                linfo("If .NO_VELOCITY. is specify, except at least entry_count = 3");
                return false;
              }
              if (!ctx.has<IOContext::NO_VELOCITY>() && entry_count < 6) {
                linfo("If .NO_VELOCITY. is NOT specify, except at least entry_count = 6");
                return false;
              }
              return true;
            }

            break;
          }
        }
      }
      current_line() = trim_spaces(m_file.get_line());
    }

    return true;
  }

  template <bool HasVelocity>
  inline bool parse_atom_line(const TokenSet& tokens, IOContext& ctx, size_t particle_count, size_t type) {

    Vec3d r{}, v{};
    size_t index = particle_count - 1;

    if constexpr (HasVelocity) {
      if (!tokens.atleast(6)) {
        linfo("ERR: line too shorts %ld %ld\n", tokens.size(), 3);
        return false;
      }
    } else {
      if (!tokens.atleast(3)) {
        linfo("ERR: line too shorts %ld %ld\n", tokens.size(), 3);
        return false;
      }
    }

    if (!parse_single_token(tokens.at<0>(), r.x))
      return false;
    if (!parse_single_token(tokens.at<1>(), r.y))
      return false;
    if (!parse_single_token(tokens.at<2>(), r.z))
      return false;

    r.x = r.x * ctx.data.cell.m11 + r.y * ctx.data.cell.m12 + r.z * ctx.data.cell.m13;
    r.x = r.x * ctx.data.cell.m21 + r.y * ctx.data.cell.m22 + r.z * ctx.data.cell.m23;
    r.x = r.x * ctx.data.cell.m31 + r.y * ctx.data.cell.m32 + r.z * ctx.data.cell.m33;

    policy::Traits::assign_position(ctx.data, index, r.x, r.y, r.z);
    policy::Traits::assign_type(ctx.data, index, type);

    if constexpr (HasVelocity) {
      if (!parse_single_token(tokens.at<0>(), v.x))
        return false;
      if (!parse_single_token(tokens.at<1>(), v.y))
        return false;
      if (!parse_single_token(tokens.at<2>(), v.z))
        return false;

      policy::Traits::assign_velocity(ctx.data, index, v.x, v.y, v.z);
    }

    return true;
  };

  template <bool HasVelocity> inline bool parse_atom_block_tmpl(IOContext& ctx) {
    double mass{};
    ssize_t type{};
    size_t particle_count{};

    policy::Traits::allocate(ctx.data, ctx.data.nat);

    while (!m_file.eof() && particle_count < ctx.data.nat) {
      if (!skipline(current_line())) {
        tokenize(current_line(), m_tokens);

        if (m_tokens.atleast(3)) {
          bool atom_overflow = cmp::gt(++particle_count, ctx.data.nat);
          if (atom_overflow || !parse_atom_line<HasVelocity>(m_tokens, ctx, particle_count, type)) {
            if (atom_overflow)
              PANIC("Too many atoms expected %ld got %ld", ctx.data.nat, particle_count);
            PANIC("Error encounterd at line: %s", std::string(current_line()).c_str())
          }
        } else {
          if (!parse_single_token(m_tokens.at<0>(), mass)) {
            ctx.species.add_species(m_tokens.at<0>());
            type = ctx.species.get_type(m_tokens.at<0>());
          }
        }
      }
      current_line() = trim_spaces(m_file.get_line());
    }

    return true;
  }

  bool parse(IOContext& ctx) override {

    current_line() = {};
    size_t entry_count{};
    if (!parse_header(ctx, entry_count))
      PANIC("Ill-formatted header !");

    current_line() = trim_spaces(m_file.get_line());
    if (ctx.has<IOContext::NO_VELOCITY>())
      parse_atom_block_tmpl<false>(ctx);
    else
      parse_atom_block_tmpl<true>(ctx);

    return true;
  }
};

} // namespace atomeye

template <> inline factory::IOFactory<Parser>::IOFactory() {
  register_format<lammps_data::LAMMPSDataParser>();
  register_format<lammps_dump::LAMMPSDumpParser>();
  register_format<xyz::ExtendedXYZParser>();
  register_format<vasp::VaspPoscarParser>();
  register_format<vasp::VaspMiladyParser>();
  register_format<atomeye::AtomEyeParser>();
}

template <> inline factory::IOFactory<Writer>::IOFactory() {
  register_format<vasp::VaspPoscarWriter>();
  register_format<vasp::VaspMiladyWriter>();
}

template <typename T> struct FileInfos {
  FileCompression compression = FileCompression::NONE;
  typename factory::IOFactory<T>::Format format;
};

template <typename T> const FileInfos<T> guess_file_infos(const std::string& filepath) {

  auto ext = std::filesystem::path(filepath).extension().string();
  if (!ext.empty() && ext.front() == '.')
    ext.erase(ext.begin()); // strip leading dot

  if (ext.empty()) {
    return {FileCompression::NONE, io_factory<T>().from_alias("")};
  }

  auto compression = convert_char_to_compression(ext);
  if (compression != FileCompression::NONE) {
    auto second = std::filesystem::path(filepath).stem().extension().string();
    if (!second.empty() && second.front() == '.')
      second.erase(second.begin()); // strip leading dot
    ext = second;
  }

  return {compression, io_factory<T>().from_alias(ext)};
}

template <typename T>
const FileInfos<T> get_file_infos(const std::string& filepath, std::string user_format, std::string user_compression) {

  auto guessed = guess_file_infos<T>(filepath);
  auto format = io_factory<T>().from_alias(user_format);

  if (!format.is_valid && !guessed.format.is_valid)
    PANIC("Could not guess file format for '%s'.", filepath.c_str());

  return {user_compression.empty() ? guessed.compression : convert_char_to_compression(user_compression),
          format.is_valid ? format : guessed.format};
}

inline bool read_external_atom_file(IOContext& ctx, const std::string& filepath, const std::string& format = "",
                                    const std::string& compression = "") {
  namespace fs = std::filesystem;
  fs::path abspath = fs::absolute(filepath);
  if (!fs::exists(abspath))
    PANIC("Input file does not exists: '%s'", filepath.c_str());

  auto infos = get_file_infos<Parser>(filepath, format, compression);
  auto parser = infos.format.creator(filepath, FileMode::READ, infos.compression);
  if (!parser)
    PANIC("Could not create atom parser for '%s'", filepath.c_str());
  bool success = (*parser)(ctx);
  return success;
}

inline bool write_external_atom_file(IOContext& ctx, const std::string& filepath, const std::string& format = "",
                                     const std::string& compression = "") {
  auto infos = get_file_infos<Writer>(filepath, format, compression);
  auto writer = infos.format.creator(filepath, FileMode::WRITE, infos.compression);
  if (!writer)
    PANIC("Could not create atom writer for '%s'", filepath.c_str());
  bool success = (*writer)(ctx);
  return success;
}

} // namespace io
} // namespace aio

/* -- OLD STUFF -- */

#ifdef OLD_STUFF
namespace policy {

struct ExaStampTraits {

  using Vec3d = typename onika::math::Vec3d;
  using Mat3d = typename onika::math::Mat3d;

  struct DataBuffer {
    size_t nat;
    size_t ntypes;
  };
};

struct InternalTraits {

  using Vec3d = typename base_types::Vec3d;
  using Mat3d = typename base_types::Mat3d;
  using IJK = typename base_types::IJK;

  struct DataBuffer {
    size_t nat;
    size_t ntypes;

    Mat3d cell{};
    Vec3d origin{};

    base_types::ParticleArray<double, 3> positions{};
    base_types::ParticleArray<double, 3> velocities{};
    base_types::ParticleArray<double, 3> forces{};
    base_types::ParticleArray<int, 1> types{};

    void reorder_atom_by_types() {
      std::vector<size_t> indices(nat);
      std::iota(indices.begin(), indices.end(), 0);

      int* __restrict types_ptr = types.data();
      std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) { return types_ptr[a] < types_ptr[b]; });

      auto permute = [&](auto& arr) {
        using ArrayT = std::decay_t<decltype(arr)>;
        using T = std::remove_reference_t<decltype(arr.data()[0])>;
        constexpr size_t DIM = ArrayT::D;
        base_types::ParticleArray<T, DIM> tmp;
        tmp.resize(nat);

        for (size_t i = 0; i < nat; ++i) {
          for (size_t d = 0; d < DIM; ++d) {
            tmp(i, d) = arr(indices[i], d);
          }
        }

        arr = std::move(tmp);
      };

      permute(positions);
    }
  };
};

#ifdef EXASTAMP_POLICY
using Traits = ExaStampTraits;
#else
using Traits = InternalTraits;
#endif

} // namespace policy

namespace matrix_traits {

using namespace base_types;

struct Mat3dSlice {
  double &x, &y, &z;
};

template <size_t I> constexpr inline double& at(Vec3d& v) {
  double* base = &v.x;
  static_assert(I >= 0 && I < 3);
  return base[I];
}

constexpr inline double& at(Vec3d& v, size_t index) {
  double* base = &v.x;
  return base[index];
}

template <size_t I, size_t J> constexpr inline double& at(Mat3d& m) {
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

inline Vec3d column(size_t i, const Mat3d& m) {
  const double* base = &m.m11;
  return Vec3d{base[i], base[i + 3], base[i + 6]};
}

inline Mat3dSlice row(size_t i, Mat3d& m) {
  double* base = &m.m11;
  size_t index = 3 * i;
  return Mat3dSlice{base[index], base[index + 1], base[index + 2]};
}

inline Vec3d row(size_t i, const Mat3d& m) {
  const double* base = &m.m11;
  size_t index = 3 * i;
  return Vec3d{base[index], base[index + 1], base[index + 2]};
}

template <typename T>
concept VecLike = requires(T v) {
  { v.x } -> std::convertible_to<double>;
  { v.y } -> std::convertible_to<double>;
  { v.z } -> std::convertible_to<double>;
};

template <VecLike V> inline void scale(V& v, double scale) {
  v.x *= scale;
  v.y *= scale;
  v.z *= scale;
}

inline void print_mat3(const Mat3d& m) {
  printf("---\n");
  for (size_t i = 0; i < 3; ++i) {
    Vec3d r = row(i, m);
    printf("%18.5e%18.5e%18.5e\n", r.x, r.y, r.z);
  }
  printf("---\n");
}

} // namespace matrix_traits

namespace io {

using namespace tokens;
using namespace tokens::numeric;
using namespace enum_traits;
using namespace basic_types;

namespace mem_size {
constexpr size_t KB = 1 << 10;
constexpr size_t MB = 1 << 20;
constexpr size_t GB = 1 << 30;

constexpr size_t READ_BUFFER_SIZE = 1 * MB;
constexpr size_t FILE_BUFFER_SIZE = READ_BUFFER_SIZE;
constexpr size_t LZMA_STREAM_BUF_SIZE = READ_BUFFER_SIZE;
constexpr size_t LZMA_INTERNAL_BUF_SIZE = READ_BUFFER_SIZE;
constexpr size_t BZ2_STREAM_BUF_SIZE = READ_BUFFER_SIZE;
constexpr size_t BZ2_INTERNAL_BUF_SIZE = READ_BUFFER_SIZE;

constexpr size_t WRITE_BUF_SIZE = 1 * MB;
constexpr size_t WRITE_BUF_THRESH = FILE_BUFFER_SIZE - 100 * KB;
} // namespace mem_size

/* ------------------------------------------------------------------------- */

struct IOContext {

  using DataBuffer = typename policy::Traits::DataBuffer;

  enum ContextFlags : uint16_t {
    REMAP_ATOM = 1 << 1,
    DOMAIN_ONLY = 1 << 2,
    TRICLINIC = 1 << 3,
  };

  uint16_t flags = 0;
  DataBuffer data{};
  Vec3d timer{};

  template <uint16_t Flag> inline bool has() { return (flags & (Flag)); }

  template <int16_t Flag> void set(bool toggle) { (toggle) ? flags |= Flag : flags &= Flag; }
};

/* ------------------------------------------------------------------------- */

inline void count_type_occurence(std::array<size_t, SpeciesMap::MAX>& counts, IOContext& ctx) {
  if (!(ctx.data.nat > 0 && ctx.data.types.size() > 0))
    return;
  int* __restrict types = ctx.data.types.data();
  for (size_t i = 0; i < ctx.data.nat; ++i) {
    counts[types[i]]++;
  }
  for (size_t i = 0; i < SpeciesMap::MAX; ++i) {
    ctx.data.species.counter.set(i, counts[i]);
  }
};

inline void count_type_occurence(IOContext& ctx) {
  std::array<size_t, SpeciesMap::MAX> counts = {0};
  count_type_occurence(counts, ctx);
};

template <typename T, size_t D>
void permute_particle_array(base_types::ParticleArray<T, D>& array, std::vector<int>& indices, size_t nat) {
  using ArrayType = std::decay_t<decltype(array)>;
  using ScalarType = std::remove_reference_t<decltype(array.data()[0])>;
  constexpr size_t dim = ArrayType::D;

  ParticleArray<T, D> tmp{};
  tmp.resize(nat);

  ScalarType* __restrict recv = tmp.data();
  const ScalarType* __restrict send = array.data();

  for (size_t i = 0; i < nat; ++i) {
    const ScalarType* __restrict src = send + indices[i] * dim;
    ScalarType* __restrict dst = recv + i * dim;
    std::memcpy(dst, src, dim * sizeof(ScalarType));
  }
  array = std::move(tmp);
}

inline void sort_atom_by_type(IOContext& ctx) {

  std::array<size_t, SpeciesMap::MAX> counts;
  std::array<size_t, SpeciesMap::MAX> offsets;
  std::vector<int> new_inds(ctx.data.nat);

  count_type_occurence(counts, ctx);

  size_t sum = 0;
  for (size_t type = 0; type < SpeciesMap::MAX; ++type) {
    offsets[type] = sum;
    sum += counts[type];
  }

  int* __restrict types = ctx.data.types.data();

  for (size_t i = 0; i < ctx.data.nat; ++i) {
    int type = types[i];
    new_inds[offsets[type]++] = i;
  }

  permute_particle_array(ctx.data.types, new_inds, ctx.data.nat);
  permute_particle_array(ctx.data.positions, new_inds, ctx.data.nat);
}

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
  virtual size_t write(const char* data, size_t size) = 0;

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
  size_t write(const char* data, size_t size) override;

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
      : m_path(std::move(path)), m_mode(mode), m_compression(compression) {}

public:
  inline const std::string& path() const { return m_path; };
  inline FileMode mode() const { return m_mode; };
  inline FileCompression compression() const { return m_compression; };
};

struct FloatFormat {
  std::chars_format fmt = std::chars_format::fixed;
  int precision = 6;
  int width = 0;
  char fill = ' ';
  bool align_right = true;
};

template <size_t BUF_SIZE>
  requires(BUF_SIZE > 1024)
class IOBuffer {
  char m_buf[BUF_SIZE];
  uint64_t m_pos = 0;

public:
  FloatFormat format{};

  static constexpr size_t threshold = BUF_SIZE - 1024;
  inline uint64_t& cursor() { return m_pos; }
  inline uint64_t size() { return m_pos; }
  inline size_t capacity() { return BUF_SIZE; }
  const char* data() const { return m_buf; }
  void clear() { m_pos = 0; }

  inline size_t cmplen(size_t a, size_t b) { return b ^ ((a ^ b) & -(a < b)); }

  template <typename T> void append(const T& value) {
    if constexpr (std::is_same_v<T, char>) {
      append_char(value);
    } else if constexpr (std::is_same_v<T, std::string_view>) {
      append_string_view(value);
    } else if constexpr (matrix_traits::VecLike<T>) {
      append_vec3d(value, format);
    } else if constexpr (std::is_floating_point_v<T>) {
      append_float(value, format);
    } else if constexpr (std::is_integral_v<T>) {
      append_int(value, format);
    } else if constexpr (std::is_convertible_v<T, std::string_view>) {
      append_string_view(std::string_view(value));
    } else {
      static_assert(std::bool_constant<false>::value, "Unsupported type");
    }
  }

  template <typename T> void append(const T& value, const FloatFormat& fmt) {
    if constexpr (std::is_floating_point_v<T>) {
      append_float(value, fmt);
    } else if constexpr (matrix_traits::VecLike<T>) {
      append_vec3d(value, fmt);
    } else if constexpr (std::is_integral_v<T>) {
      append_int(value, fmt);
    } else {
      static_assert(std::bool_constant<false>::value, "Unsupported type");
    }
  };

  void append_char(char value) {
    if (m_pos >= BUF_SIZE)
      return;
    m_buf[m_pos++] = value;
  }

  void append_string_view(StringView value) {
    size_t len = cmplen(value.size(), BUF_SIZE - m_pos);
    std::memcpy(m_buf + m_pos, value.data(), len);
    m_pos += len;
  }

  template <typename T>
    requires std::is_floating_point_v<T>
  inline void append_float(T value, const FloatFormat& fmt) {
    char* begin = m_buf + m_pos;
    char* end = begin + 32;
    auto [ptr, ec] = std::to_chars(begin, end, value, fmt.fmt, fmt.precision);
    if (ec != std::errc()) [[unlikely]]
      return;

    size_t len = static_cast<size_t>(ptr - begin);

    if (fmt.width > 0 && len < static_cast<size_t>(fmt.width)) {
      size_t pad = fmt.width - len;
      if (fmt.align_right) {
        std::memmove(begin + pad, begin, len);
        std::fill_n(begin, pad, fmt.fill);
      } else {
        std::fill_n(ptr, pad, fmt.fill);
      }
      len = fmt.width;
    }
    m_pos += len;
  }

  template <typename T>
    requires std::is_integral_v<T>
  inline void append_int(T value, const FloatFormat& fmt) {
    char* begin = m_buf + m_pos;
    char* end = begin + 32;
    auto [ptr, ec] = std::to_chars(begin, end, value);
    if (ec != std::errc()) [[unlikely]]
      return;

    size_t len = static_cast<size_t>(ptr - begin);

    if (fmt.width > 0 && len < static_cast<size_t>(fmt.width)) {
      size_t pad = fmt.width - len;
      if (fmt.align_right) {
        std::memmove(begin + pad, begin, len);
        std::fill_n(begin, pad, fmt.fill);
      } else {
        std::fill_n(ptr, pad, fmt.fill);
      }
      len = fmt.width;
    }
    m_pos += len;
  }

  template <matrix_traits::VecLike V> inline void append_vec3d(V value, const FloatFormat& fmt) {
    append(value.x, fmt);
    space();
    append(value.y, fmt);
    space();
    append(value.z, fmt);
  };

  inline constexpr void space() { append(' '); }
  inline constexpr void newline() { append('\n'); }
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

  template <size_t N> inline void write(IOBuffer<N>& buf, bool force = false) {
    if (m_mode != FileMode::WRITE)
      PANIC("Can't write, the file is not opened in write mode.");
    if (buf.size() > IOBuffer<N>::threshold || force) {
      m_handler->write(buf.data(), buf.size());
      buf.clear();
    }
  }
};

/* ------------------------------------------------------------------------- */

struct Metadata {
  const std::string name = "";
  const std::vector<std::string> aliases = {};
  const std::string description = "";
};

template <typename Derived> struct IOTraitsBase {
  static Metadata metadata() {
    return {.name = std::string(Derived::name),
            .aliases = std::vector<std::string>(Derived::aliases.begin(), Derived::aliases.end()),
            .description = std::string(Derived::description)};
  }
};

template <typename T>
concept HasValidMetadataTraits = requires {
  typename T::Traits;
  { T::Traits::name } -> std::convertible_to<const char*>;
  {
    T::Traits::aliases
  } -> std::convertible_to<const std::array<const char*, std::tuple_size<decltype(T::Traits::aliases)>::value>&>;
  { T::Traits::description } -> std::convertible_to<const char*>;
  { T::Traits::metadata() } -> std::same_as<Metadata>;

  requires T::Traits::name[0] != '\0';
  requires T::Traits::description[0] != '\0';
  requires T::Traits::aliases[0][0] != '\0';
};

template <typename T>
  requires HasValidMetadataTraits<T>
const Metadata GetMetadata() {
  return T::Traits::metadata();
}

/* ------------------------------------------------------------------------- */

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

class Writer {
public:
  Writer() = default;
  virtual ~Writer() = default;
  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;
  Writer(Writer&&) = delete;
  Writer& operator=(Writer&&) = delete;

  virtual bool write(IOContext& ctx) = 0;
  virtual bool operator()(io::IOContext&) { return false; };
};

class TextWriter : public Writer {
protected:
  TextFile m_file;
  IOBuffer<mem_size::WRITE_BUF_SIZE> m_buf;

public:
  TextWriter(std::string filepath, FileMode mode, FileCompression compression)
      : m_file(std::move(filepath), mode, compression) {}
  virtual ~TextWriter() override = default;

  // bool write(IOContext&) override { std::cout << "ok" << std::endl; return false; };
  bool operator()(IOContext& ctx) override;
};

/* ------------------------------------------------------------------------- */

template <typename Base> class IOFactory {
public:
  using creator_t = std::unique_ptr<Base> (*)(std::string, FileMode, FileCompression);

  struct Format {
    const Metadata metadata;
    const creator_t creator;
    bool is_valid = true;
  };

  IOFactory();

  std::unordered_map<std::string, size_t> alias_to_index;
  std::unordered_map<size_t, Format> index_to_parser;

  template <typename T>
  static std::unique_ptr<Base> create_instance(std::string filepath, FileMode mode, FileCompression compression) {
    return std::make_unique<T>(std::move(filepath), mode, compression);
  }

  template <typename T> void register_format() {

    size_t index = index_to_parser.size();
    Metadata metadata = GetMetadata<T>();
    index_to_parser.insert({index, {.metadata = metadata, .creator = &create_instance<T>}});

    // map extension to the index
    for (const std::string& alias : metadata.aliases) {
      std::unordered_map<std::string, size_t>::const_iterator it = alias_to_index.find(alias);
      if (it != alias_to_index.end()) {
        continue; // the alias is already used
      }
      alias_to_index.insert({alias, index});
    }
  }

  Format from_alias(const std::string& alias) const {
    std::unordered_map<std::string, size_t>::const_iterator it = alias_to_index.find(alias);
    if (it == alias_to_index.end()) {
      return {.metadata = Metadata{}, .creator = nullptr, .is_valid = false};
    }
    return index_to_parser.find(it->second)->second;
  }

  inline std::string infos() const {
    std::string infos = "Supported atom format:\n";
    // typename decltype(index_to_parser)::const_iterator it;
    size_t size = index_to_parser.size();
    for (size_t i = 0; i < size; ++i) {
      infos += strf("\t- %s\n", index_to_parser.find(i)->second.metadata.name.c_str());
    }
    infos += "\n";
    return infos;
  }
};

template <typename T> const IOFactory<T>& io_factory() {
  static const IOFactory<T> instance{};
  return instance;
}

/* ------------------------------------------------------------------------- */
// Parser for LAMMPS Data format

namespace lammps_data {
using namespace tokens;
using namespace tokens::numeric;

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

constexpr std::array<size_t, 1> SIZE_KEYWORD_ATOMS{NEEDLES_KEYWORD_ATOMS[0].size()};
constexpr std::array<size_t, 2> SIZE_KEYWORD_ATOM_TYPES{NEEDLES_KEYWORD_ATOM_TYPES[0].size(),
                                                        NEEDLES_KEYWORD_ATOM_TYPES[1].size()};
constexpr std::array<size_t, 2> SIZE_KEYWORD_XAXIS{NEEDLES_KEYWORD_XAXIS[0].size(), NEEDLES_KEYWORD_XAXIS[1].size()};
constexpr std::array<size_t, 2> SIZE_KEYWORD_YAXIS{NEEDLES_KEYWORD_YAXIS[0].size(), NEEDLES_KEYWORD_YAXIS[1].size()};
constexpr std::array<size_t, 2> SIZE_KEYWORD_ZAXIS{NEEDLES_KEYWORD_ZAXIS[0].size(), NEEDLES_KEYWORD_ZAXIS[1].size()};

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
    "Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals",
};

bool parse_keyword_atoms(const TokenSet&, IOContext&);
bool parse_keyword_atom_types(const TokenSet&, IOContext&);
bool parse_keyword_xaxis(const TokenSet&, IOContext&);
bool parse_keyword_yaxis(const TokenSet&, IOContext&);
bool parse_keyword_zaxis(const TokenSet&, IOContext&);
bool parse_keyword_tilts(const TokenSet&, IOContext&);
bool parse_keyword_section(const TokenSet&, IOContext&);

constexpr std::array<KeywordParserFn, enum_size<Keyword>()> keyword_parsers{
    &parse_keyword_atoms, &parse_keyword_atom_types, &parse_keyword_xaxis,   &parse_keyword_yaxis,
    &parse_keyword_zaxis, &parse_keyword_tilts,      &parse_keyword_section,
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

  // ParticleTupleIO& p = ctx.data.particles[index];
  // p[field::rx] = x;
  // p[field::ry] = y;
  // p[field::rz] = z;
  // p[field::vx] = 0.0;
  // p[field::vy] = 0.0;
  // p[field::vz] = 0.0;
  // p[field::id] = index;
  // p[field::type] = type - 1;

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

  // ParticleTupleIO& p = ctx.data.particles[index];
  // p[field::rx] = x;
  // p[field::ry] = y;
  // p[field::rz] = z;
  // p[field::vx] = 0.0;
  // p[field::vy] = 0.0;
  // p[field::vz] = 0.0;
  // p[field::id] = index;
  // p[field::type] = type - 1;

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

  // ParticleTupleIO& p = ctx.data.particles[index];
  // p[field::vx] = vx;
  // p[field::vy] = vy;
  // p[field::vz] = vz;
  return true;
}

struct LAMMPSDataTraits : IOTraitsBase<LAMMPSDataTraits> {
  static constexpr const char* name = "LAMMPS Data";
  static constexpr std::array<const char*, 3> aliases{"data", "lmp-data", "lmp"};
  static constexpr const char* description = "LAMMPS Data format";
};

class LAMMPSDataParser : public TextParser {
public:
  using Traits = LAMMPSDataTraits;
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
    "id", "mol", "type", "element", "x",   "y",   "z",  "xu", "yu", "zu",
    "xs", "ys",  "zs",   "xsu",     "ysu", "zsu", "vx", "vy", "vz",
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
    v.z = 0.0;
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

  ctx.data.cell.m21 = 0.0;
  ctx.data.cell.m31 = 0.0;
  ctx.data.cell.m32 = 0.0;

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

template <bool HasMolID> inline bool parse_field_molid(size_t& molid, const TokenSet& tokens, size_t index) {
  if constexpr (HasMolID) {
    return parse_single_token(tokens[index], molid);
  } else {
    molid = 0;
    return true;
  }
}

template <bool HasParticleType> inline bool parse_field_type(size_t& type, const TokenSet& tokens, size_t index) {
  if constexpr (HasParticleType) {
    return parse_single_token(tokens[index], type);
  } else {
    type = 0;
    return true;
  }
}

template <bool HasVelocities> inline bool parse_field_velocities(Vec3d& v, const TokenSet& tokens, const IJK& inds) {
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

  double* positions = ctx.data.positions.data();
  positions[index * 3 + 0] = pos.x;
  positions[index * 3 + 1] = pos.y;
  positions[index * 3 + 2] = pos.z;

  double* velocities = ctx.data.velocities.data();
  velocities[index * 3 + 0] = vel.x;
  velocities[index * 3 + 1] = vel.y;
  velocities[index * 3 + 2] = vel.z;

  ctx.data.types.data()[index] = static_cast<int>(type - 1);

  // ParticleTupleIO& p = ctx.data.particles[index];
  // p[field::rx] = pos.x;
  // p[field::ry] = pos.y;
  // p[field::rz] = pos.z;
  // p[field::vx] = vel.x;
  // p[field::vy] = vel.y;
  // p[field::vz] = vel.z;
  // p[field::id] = index;
  // p[field::type] = type - 1;

  return true;
}

struct LAMMPSDumpTraits : IOTraitsBase<LAMMPSDumpTraits> {
  static constexpr const char* name = "LAMMPS Dump";
  static constexpr std::array<const char*, 3> aliases{"dump", "lmp-dump", "lammps-dump"};
  static constexpr const char* description = "LAMMPS Dump format";
};

class LAMMPSDumpParser : public TextParser {
public:
  using Traits = LAMMPSDumpTraits;
  LAMMPSDumpParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override;
  bool parse(IOContext& ctx) override;
};

} // namespace lammps_dump

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
    "species", "pos", "velo", "id", "type",
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

template <bool HasMolID> inline bool parse_field_molid(size_t& molid, const TokenSet& tokens, size_t index) {
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
    type = ctx.data.species.get(tokens[i]);
    return cmp::ne(type, -1);
  } else if constexpr (HasParticleType) {
    return parse_single_token(tokens[j], type);
  } else {
    type = 0;
    return true;
  }
}

template <bool HasVelocities> inline bool parse_field_velocities(Vec3d& v, const TokenSet& tokens, const IJK& inds) {
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
  size_t ncol = 0;
  size_t ifield = 0;
  size_t itype = 0;

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

  // ParticleTupleIO& p = ctx.data.particles[index];
  // p[field::rx] = pos.x;
  // p[field::ry] = pos.y;
  // p[field::rz] = pos.z;
  // p[field::vx] = vel.x;
  // p[field::vy] = vel.y;
  // p[field::vz] = vel.z;
  // p[field::id] = index;
  // p[field::type] = static_cast<size_t>(type - 1);

  return true;
}

struct ExtendedXYZTraits : IOTraitsBase<ExtendedXYZTraits> {
  static constexpr const char* name = "Extended XYZ";
  static constexpr std::array<const char*, 3> aliases{"xyz", "exyz", "ext-xyz"};
  static constexpr const char* description = "Extended XYZ";
};

class ExtendedXYZParser : public TextParser {
public:
  using Traits = ExtendedXYZTraits;
  ExtendedXYZParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}

  int64_t next() override;
  bool parse(IOContext& ctx) override;

private:
  TokenSet m_ptokens;
  Properties m_properties;
};

} // namespace xyz

/* ------------------------------------------------------------------------- */

namespace vasp {

constexpr TokenNeedles<1> CARTESIAN_COORDINATES{"Cartesian"};

template <bool Cartesian>
inline bool parse_poscar_atom_line(const TokenSet& tokens, IOContext& ctx, size_t particle_count, const Vec3d& scale) {

  size_t index = particle_count - 1;
  Vec3d pos;

  if (!tokens.atleast(3))
    return false;
  if (!parse_single_token(tokens.at<0>(), pos.x))
    return false;
  if (!parse_single_token(tokens.at<1>(), pos.y))
    return false;
  if (!parse_single_token(tokens.at<2>(), pos.z))
    return false;

  double* positions = ctx.data.positions.data();

  if constexpr (Cartesian) {
    positions[index * 3 + 0] = pos.x * scale.x;
    positions[index * 3 + 1] = pos.y * scale.y;
    positions[index * 3 + 2] = pos.z * scale.z;
  } else {
    positions[index * 3 + 0] = pos.x * ctx.data.cell.m11 + pos.y * ctx.data.cell.m12 + pos.z * ctx.data.cell.m13;
    positions[index * 3 + 1] = pos.x * ctx.data.cell.m21 + pos.y * ctx.data.cell.m22 + pos.z * ctx.data.cell.m23;
    positions[index * 3 + 2] = pos.x * ctx.data.cell.m31 + pos.y * ctx.data.cell.m32 + pos.z * ctx.data.cell.m33;
  }

  ctx.data.types.data()[index] = ctx.data.species.counter.from_index(index);
  return true;
}

template <bool Milady> struct VaspTraits : IOTraitsBase<VaspTraits<Milady>> {};

template <> struct VaspTraits<false> : IOTraitsBase<VaspTraits<false>> {
  static constexpr const char* name = "VASP POSCAR";
  static constexpr std::array<const char*, 3> aliases{"vasp", "poscar", "pos"};
  static constexpr const char* description = "???";
};

template <> struct VaspTraits<true> : IOTraitsBase<VaspTraits<true>> {
  static constexpr const char* name = "VASP POSCAR MILADY";
  static constexpr std::array<const char*, 2> aliases{"mld", "milady"};
  static constexpr const char* description = "???";
};

template <bool Milady> class VaspPoscarParser : public TextParser {
public:
  using Traits = VaspTraits<Milady>;
  VaspPoscarParser(std::string filepath, FileMode mode, FileCompression compression)
      : TextParser(filepath, mode, compression) {}
  int64_t next() override;
  bool parse(IOContext& ctx) override;
};

template <bool Milady> class VaspPoscarWriter : public TextWriter {
public:
  using Traits = VaspTraits<Milady>;
  VaspPoscarWriter(std::string filepath, FileMode mode, FileCompression compression)
      : TextWriter(filepath, mode, compression) {}
  bool write(IOContext& ctx) override {

    // Vasp requires same type to be contiguous
    if (ctx.data.species.count() > 1) {
      sort_atom_by_type(ctx);
    } else {
      count_type_occurence(ctx);
    }

    IOBuffer<mem_size::WRITE_BUF_SIZE> buf;
    buf.format = {.width = 15, .fill = ' '};

    // Comment line
    buf.append("111 ");
    buf.append(ctx.data.species.count(), {.width = 0});
    for (size_t i = 0; i < ctx.data.species.count(); ++i) {
      buf.space();
      buf.append(ctx.data.species.get(i));
      buf.space();
      buf.append(0, {.width = 0});
    }
    buf.append(" 0.00 0.00 0.00\n");

    // Scaling Factor
    buf.append(1.0);
    buf.newline();
    // Lattice vectors
    buf.append(matrix_traits::column(0, ctx.data.cell));
    buf.newline();
    buf.append(matrix_traits::column(1, ctx.data.cell));
    buf.newline();
    buf.append(matrix_traits::column(2, ctx.data.cell));
    buf.newline();
    // Species
    for (size_t i = 0; i < ctx.data.species.count(); ++i) {
      buf.space();
      buf.append(ctx.data.species.get(i));
    }
    buf.newline();
    // Species count
    for (size_t i = 0; i < ctx.data.species.count(); ++i) {
      buf.space();
      buf.append(ctx.data.species.counter.get(i), {});
    }
    buf.newline();
    // write header to buffer
    m_file.write(buf, true);

    // Positions
    buf.append("Cartesian\n");
    double* __restrict positions = ctx.data.positions.data();
    for (size_t i = 0; i < ctx.data.nat; ++i) {
      buf.append(positions[i * 3 + 0]);
      buf.space();
      buf.append(positions[i * 3 + 1]);
      buf.space();
      buf.append(positions[i * 3 + 2]);
      buf.newline();
      m_file.write(buf);
    }

    if constexpr (Milady) {
      // Forces
      buf.newline();

      if (ctx.data.forces.size() == ctx.data.nat && false) {
        double* __restrict forces = ctx.data.forces.data();
        for (size_t i = 0; i < ctx.data.nat; ++i) {
          buf.append(forces[i * 3 + 0]);
          buf.space();
          buf.append(forces[i * 3 + 1]);
          buf.space();
          buf.append(forces[i * 3 + 2]);
          buf.newline();
          m_file.write(buf);
        }
      } else {
        for (size_t i = 0; i < ctx.data.nat; ++i) {
          buf.append(0.0);
          buf.space();
          buf.append(0.0);
          buf.space();
          buf.append(0.0);
          buf.newline();
          m_file.write(buf);
        }
      }

      // Stress and spin flag
      buf.newline();
      buf.append("0.0 0.0 0.0 0.0 0.0 0.0\n\n0\n");
    }

    m_file.write(buf, true);
    return true;
  }
};

} // namespace vasp

/* ------------------------------------------------------------------------- */

std::string remove_leading_chars(const std::string& str, const char c);

template <typename T> struct FileInfos {
  FileCompression compression = FileCompression::NONE;
  const typename IOFactory<T>::Format format;
};

template <typename T> const FileInfos<T> guess_file_infos(const std::string& filepath) {

  std::string extension;
  std::filesystem::path path(filepath);

  // remove leading dots return by extension()
  extension = remove_leading_chars(path.extension(), '.');

  // if no extension can't do anything here
  if (extension.empty()) {
    return {.compression = FileCompression::NONE, .format = io_factory<T>().from_alias("")};
  }

  FileCompression compression = convert_char_to_compression(extension);

  // if first extension is related to compression, get the second one
  if (compression != FileCompression::NONE) {
    extension = remove_leading_chars(path.stem().extension(), '.');
  }

  // return file infos either if format is invalid, this is handled after
  return {.compression = compression, .format = io_factory<T>().from_alias(extension)};
}

template <typename T>
const FileInfos<T> get_file_infos(const std::string& filepath, std::string user_format, std::string user_compression) {

  // never trust the user input
  FileInfos<T> guessed_infos = guess_file_infos<T>(filepath);
  const typename IOFactory<T>::Format format = io_factory<T>().from_alias(user_format);

  if ((!format.is_valid) && (!guessed_infos.format.is_valid)) {
    PANIC("Could not guess file format.");
  }

  FileInfos<T> infos{
      .compression =
          user_compression.empty() ? guessed_infos.compression : convert_char_to_compression(user_compression),
      .format = format.is_valid ? format : guessed_infos.format,
  };

  return infos;
}

bool read_external_atom_file(IOContext& ctx, std::string filepath, std::string format, std::string compression);
bool dump_external_atom_file(IOContext& ctx, std::string filepath, std::string format, std::string compression);

} // namespace io
#endif
