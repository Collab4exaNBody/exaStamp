#pragma once

#include <cctype>
#include <charconv>
#include <cstring>
#include <span>
#include <string_view>

namespace tokens {

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

struct ColumnDelimiter {
  constexpr bool operator()(char c) const noexcept { return c == ':'; }
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

// template <typename Delimiter>
// inline void tokenizer_tmpl(const std::string_view str, TokenSet& tokens, Delimiter&& delimiter) {
//   tokens.clear(); // reset the token set
//   split(str, delimiter, [&tokens](std::string_view token) {
//     if (!token.empty() && token.front() != '#') {
//       tokens.push_back(token);
//     }
//   });
// }

template <typename Delimiter, size_t N>
inline void tokenizer_tmpl(const std::string_view str, TokenSetTmpl<N>& tokens, Delimiter&& delimiter) {
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

inline bool starts_with(std::string_view line, std::string_view prefix) {
  return ((line.size() >= prefix.size()) && (std::memcmp(line.data(), prefix.data(), prefix.size()) == 0));
}

}
