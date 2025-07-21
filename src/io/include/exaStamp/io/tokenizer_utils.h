#pragma once

#include <charconv>
#include <cmath>
#include <cstdarg>
#include <cstring>
#include <span>
#include <string>
#include <string_view>

#define DEFAULT_TOKEN_SET_SIZE 128
#define SPLIT_MAX_COUNT 1000

inline std::string strf(const char* format, ...) {
  char buffer[1024];
  va_list args;
  va_start(args, format);
  vsprintf(buffer, format, args);
  return std::string(buffer);
};

template <typename T, typename C> inline constexpr T enum_to(C EnumMember) { return static_cast<T>(EnumMember); }

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

namespace tokens {

using StringView = std::string_view;
using StringViewSpan = std::span<const StringView>;

inline bool startswith(StringView str, StringView prefix) {
  return (str.size() >= prefix.size()) && (std::memcmp(str.data(), prefix.data(), prefix.size()) == 0);
}

inline bool endswith(StringView str, StringView suffix) {
  return (str.size() >= suffix.size()) &&
         (std::memcmp(str.data() + (str.size() - suffix.size()), suffix.data(), suffix.size()) == 0);
}

/* ------------------------------------------------------------------------- */

// Struct that holds a fixed size array of token.
// Warning: Tokens are std::string_view so the TokenSet does not own the data
template <size_t Size> struct TokenSetTmpl {
  alignas(64) std::array<StringView, Size> tokens{};
  size_t len = 0;

  inline void clear() noexcept { len = 0; }
  inline void reset() noexcept {
    tokens = {0};
    clear();
  }

  inline StringView* data() noexcept { return tokens.data(); }
  inline const StringView* data() const noexcept { return tokens.data(); }

  inline StringView& operator[](size_t i) noexcept { return tokens[i]; }
  inline const StringView& operator[](size_t i) const noexcept { return tokens[i]; }

  template <size_t index> inline constexpr StringView& at() noexcept {
    static_assert(index < Size);
    return tokens[index];
  }

  template <size_t index> inline constexpr const StringView& at() const noexcept {
    static_assert(index < Size);
    return tokens[index];
  }

  constexpr size_t size() const noexcept { return len; }
  constexpr size_t capacity() const noexcept { return Size; }
  constexpr inline bool full() const noexcept { return size() >= capacity(); }
  constexpr inline bool atleast(size_t n) const noexcept { return size() >= n; }

  inline void push_back(StringView token) noexcept {
    if (!full()) {
      tokens[len++] = token;
    }
  }

  inline void push_back_unchecked(StringView token) noexcept { tokens[len++] = token; }

  // iteration accessor
  StringView* begin() noexcept { return tokens.begin(); }
  const StringView* begin() const noexcept { return tokens.begin(); }
  StringView* end() noexcept { return tokens.begin() + len; }
  const StringView* end() const noexcept { return tokens.begin() + len; }
};

using TokenSet = TokenSetTmpl<DEFAULT_TOKEN_SET_SIZE>;

/* ------------------------------------------------------------------------- */

template <char... Chars> struct IsChar {
  constexpr bool operator()(char c) const noexcept { return ((c == Chars) || ...); }
  constexpr bool operator()(const char* cptr) const noexcept { return cptr && ((*cptr == Chars) || ...); }
  // constexpr bool operator()(const char* cptr) const noexcept { return cptr && this(*cptr); } //
  // Equivalent ?
};

struct IsSpace : IsChar<' ', '\t', '\n', '\r'> {};
struct SpaceDelimiter : IsChar<' ', '\t'> {};
struct CommaDelimiter : IsChar<','> {};
struct ColumnDelimiter : IsChar<':'> {};
struct SemiColumnDelimiter : IsChar<';'> {};
struct DoubleQuoteDelimiter : IsChar<'"'> {};
struct SingleQuoteDelimiter : IsChar<'\''> {};
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
  while (!str.empty() && predicate(str.front())) {
    str.remove_prefix(1);
  }
  while (!str.empty() && predicate(str.back())) {
    str.remove_suffix(1);
  }
  return str;
}

inline StringView trim_spaces(StringView str) noexcept { return trim(str, IsSpace{}); }

template <typename Predicate> inline bool skipline_tmpl(StringView str, Predicate&& predicate) noexcept {
  return str.empty() || predicate(str.front());
}

inline bool skipline(StringView str) noexcept { return skipline_tmpl(str, IsHashComment{}); }

template <typename Predicate> inline bool front(StringView str, Predicate&& predicate) { return predicate(str); };

template <char... Chars> inline bool front(StringView str) { return front(str, IsChar<Chars...>{}); }

template <typename Callback, typename Delimiter>
inline void split(StringView str, Delimiter&& delimiter, Callback&& callback,
                  size_t max_count = SPLIT_MAX_COUNT) noexcept {
  size_t count = 0;
  const char* data = str.data();
  const char* end = data + str.size();
  while (data < end && cmp::lt(count, max_count)) {
    while (data < end && delimiter(data)) {
      ++data;
    }
    const char* token_start = data;
    while (data < end && !delimiter(data)) {
      ++data;
    }
    if (token_start != data) {
      callback(StringView(token_start, static_cast<size_t>(data - token_start)));
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

// Return true if 'token' match any of the needles
// Stop at the first match and set 'index' to the matching index in 'needles'.
inline bool token_match_any(StringView token, StringViewSpan needles, size_t& index) {
  const size_t nc = needles.size();
  if (nc == 0) {
    return false;
  }
  for (size_t i = 0; i < nc; ++i) {
    if (token == needles[i]) {
      index = i;
      return true;
    }
  }
  return false;
}

inline bool token_match_any(StringView token, StringViewSpan needles) {
  size_t index;
  return token_match_any(token, needles, index);
}

// Return true when any token in 'tokens' match any needle in 'needles'
// Stop at the first match and provide i and j indices of the token and needles respectively.
inline bool token_set_match_any(StringViewSpan tokens, StringViewSpan needles, size_t& ii, size_t& jj) {
  const size_t token_count = tokens.size();
  const size_t needle_count = needles.size();

  if (needle_count == 0) {
    return false;
  }

  for (size_t i = 0; i < token_count; ++i) {
    for (size_t j = 0; j < needle_count; ++j) {
      if (tokens[i] == needles[j]) {
        ii = i; // token index
        jj = j; // needle index
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
inline bool token_set_match_seq(StringViewSpan tokens, StringViewSpan needles, size_t& index) {
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

template <size_t N> using TokenNeedles = std::array<StringView, N>;

template <size_t N> inline bool token_match_any(const StringView token, const TokenNeedles<N>& needles) {
  return token_match_any(token, StringViewSpan(needles.data(), N));
}

template <size_t N> inline bool token_match_any(const StringView token, const TokenNeedles<N>& needles, size_t& index) {
  return token_match_any(token, StringViewSpan(needles.data(), N), index);
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

/* ------------------------------------------------------------------------- */
// convertion functions

namespace numeric {

// template <typename T> bool parse_single_token_impl(const StringView token, T& rval) {
//   static_assert(std::is_arithmetic_v<T>, "parse_single_token_impl requires arithmetic type");
//   auto [ptr, err] = std::from_chars(token.begin(), token.begin() + token.size(), static_cast<T&>(rval));
//   return (err == std::errc() && ptr == token.end());
// }

template <typename T, typename = void>
struct TokenParser;

template <typename T>
struct TokenParser<T, std::enable_if_t<std::is_floating_point_v<T>>> {
  static bool parse(const char* begin, size_t size, T& rval) {
    const char* end = begin + size;
    auto [ptr, err] = std::from_chars(begin, end, static_cast<T&>(rval), std::chars_format::general);
    return (err == std::errc() && ptr == end);
  }
};

template <typename T>
struct TokenParser<T, std::enable_if_t<std::is_integral_v<T>>> {

  static constexpr T min = std::numeric_limits<T>::min();
  static constexpr T max = std::numeric_limits<T>::max();

  static bool parse(const char* begin, size_t size, T& rval) {
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

template <typename T>
inline bool parse_single_token(const StringView token, T& rval) {
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

/* ----------------------- Legacy version ----------------------- */
// convert a string_view token to numeric type
// template <typename T> bool token_to_num(const StringView token, T& rval) {
//   auto [ptr, err] = std::from_chars(token.begin(), token.begin() + token.size(), rval);
//   return (err == std::errc() && ptr == token.end());
// }

/***
// Dark magic that expand to tokens_to_num with automatic type infering
#define parb ()
#define IOexpand(...) IOexpandB(IOexpandB(IOexpandB(IOexpandB(__VA_ARGS__))))
#define IOexpandB(...) IOexpandA(IOexpandA(IOexpandA(IOexpandA(__VA_ARGS__))))
#define IOexpandA(...) __VA_ARGS__
#define call_parse(token, value) token_to_num<decltype(value)>(token, value)
#define IOfor_each(macro, ...) __VA_OPT__(IOexpand(IOfor_each_helper(macro, __VA_ARGS__)))
#define IOfor_each_helper(macro, token, value, ...)                                                \
  macro(token, value) __VA_OPT__(IOfor_each_again parb(macro, __VA_ARGS__))
#define IOfor_each_again() &&IOfor_each_helper
#define tokens_to_num(...) IOfor_each(call_parse, __VA_ARGS__)
***/
/* -----------------------  Legacy version ----------------------- */

} // namespace tokens

namespace enum_traits {

template <typename T>
concept has_count_member = requires { T::Count; };

template <typename T>
concept HasMemberCount = requires { T::Count; };

template <typename T>
concept HasMemberNone = requires { T::None; };

template <has_count_member T> constexpr size_t enum_size() { return static_cast<size_t>(T::Count); }

template <typename T> constexpr size_t enum_to_index(T t) { return static_cast<size_t>(t); }

template <typename T> constexpr T enum_from_index(size_t index) { return static_cast<T>(index); }

template <typename Enum, typename Type> struct BasicEnumContainer {

  constexpr inline void set(Enum e, const Type& value) {
    static_assert(std::is_enum<Enum>());
    if (enum_to_index(e) > enum_size<Enum>()) {
      return;
    }
    _data[enum_to_index(e)] = value;
  }

  template <Enum e> constexpr inline void set(const Type& value) { set(e, value); }

  inline const Type& operator[](size_t index) const noexcept { return _data[index]; }
  inline const Type& operator[](Enum e) const noexcept { return _data[enum_to_index(e)]; }

private:
  std::array<Type, enum_size<Enum>()> _data{};
};

template <typename Enum, typename Type, typename... Args> struct Builder;

template <typename Enum, typename Type, typename Value, typename... Args>
  requires std::is_convertible_v<Value, Type>
struct Builder<Enum, Type, Enum, Value, Args...> {
  static constexpr void build(BasicEnumContainer<Enum, Type>& container, Enum e, Value v, Args... rest) {
    container.set(e, Type(v));
    Builder<Enum, Type, Args...>::build(container, rest...);
  }
};

template <typename Enum, typename Type> struct Builder<Enum, Type> {
  static constexpr void build(BasicEnumContainer<Enum, Type>&) {}
};

template <typename Enum, typename Type, typename... Args>
constexpr BasicEnumContainer<Enum, Type> make_enum_container(Args... args) {
  static_assert(sizeof...(args) % 2 == 0, "Must provide pairs of (Enum, Type)");
  BasicEnumContainer<Enum, Type> container{};
  Builder<Enum, Type, Args...>::build(container, args...);
  return container;
}

} // namespace enum_traits

namespace properties_traits {
using namespace enum_traits;

template <typename Enum>
concept IsEnum = std::is_enum_v<Enum>;

template <typename Type, typename Enum, Enum Type::* Member>
concept TypeHasEnumMember = requires(Type t) {
  { t.*Member } -> std::same_as<Enum&>;
};

template <typename Type, IsEnum Enum, Enum Type::* Member>
  requires TypeHasEnumMember<Type, Enum, Member> && HasMemberCount<Enum> && HasMemberNone<Enum>
struct GenericContainer {

  using enum_t = Enum;

protected:
  size_t length = 0;
  std::array<Type, enum_size<Enum>()> _data{};

public:
  template <Enum e> constexpr const Type& get() const {
    static_assert(enum_to_index(e) < enum_size<Enum>(), "PropertiesContainer::get : out-of-bound");
    return _data[enum_to_index(e)];
  }

  template <Enum... args> constexpr bool has() const { return ((get<args>().*Member == args) && ...); }

  template <Enum... args> constexpr bool any() const { return ((get<args>().*Member == args) || ...); }

  template <typename... Args> constexpr void expect(Args... args) {
    static_assert((std::is_same_v<Args, Enum> && ...), "Only E enum type allowed");
    ((length += (has(args)) ? static_cast<size_t>(1) : static_cast<size_t>(0)), ...);
  }

  template <Enum... args> constexpr void expect() {
    // static_assert((std::is_same_v<Enum, Enum> && args...), "Only E enum type allowed");
    ((length += (has<args>()) ? static_cast<size_t>(1) : static_cast<size_t>(0)), ...);
  }

  const Type& get(Enum e) const {
    if (enum_to_index(e) < enum_size<Enum>()) {
      return _data[enum_to_index(e)];
    } else {
      return _data[enum_to_index(Enum::None)];
    }
  }

  template <typename... Args> bool has(Args... args) {
    static_assert((std::is_same_v<Args, Enum> && ...), "Only E enum type allowed");
    return ((get(args).*Member == args) && ...);
  }

  void set(Type t) {
    if (enum_to_index<Enum>(t.*Member) < enum_size<Enum>()) {
      _data[enum_to_index<Enum>(t.*Member)] = t;
    }
  }

  inline void clear() {
    length = 0;
    std::fill(_data.begin(), _data.end(), Type{});
  }

  inline size_t len(void) const { return length; }
  inline size_t len(void) { return length; }
};
} // namespace properties_traits

namespace functor {

template <template <bool...> typename FunctorType, size_t N> struct BooleanDispatcher {
private:
  template <size_t I, typename Seq> struct InstanceHolder;

  template <size_t I, size_t... Idx> struct InstanceHolder<I, std::index_sequence<Idx...>> {

    static constexpr std::array<bool, N> flags = [] {
      std::array<bool, N> bits{};
      for (size_t j = 0; j < N; ++j) {
        bits[N - 1 - j] = (I >> j) & 1;
      }
      return bits;
    }();

    static constexpr FunctorType<flags[Idx]...> instance{};
  };

  template <typename IndexSequence> struct FunctorTypeHelper;

  template <size_t... Idx> struct FunctorTypeHelper<std::index_sequence<Idx...>> {
    using type = FunctorType<(static_cast<void>(Idx), true)...>;
  };

  using FunctorBasePtr = const typename FunctorTypeHelper<std::make_index_sequence<N>>::type::Base*;

  template <size_t I> static constexpr FunctorBasePtr make_instance() {
    return &InstanceHolder<I, std::make_index_sequence<N>>::instance;
  }

  template <size_t... Is> static constexpr auto generate_table_impl(std::index_sequence<Is...>) {
    return std::array<FunctorBasePtr, sizeof...(Is)>{make_instance<Is>()...};
  }

  static constexpr auto generate_table() { return generate_table_impl(std::make_index_sequence<(1 << N)>{}); }

  inline static constexpr std::array<FunctorBasePtr, (1 << N)> table = generate_table();

public:
  template <typename... Args> static const FunctorBasePtr dispatch(Args... args) {
    static_assert(sizeof...(Args) == N);
    static_assert((std::is_constructible_v<Args, bool> && ...));
    size_t index = 0;
    bool bits[] = {static_cast<bool>(args)...};
    for (size_t i = 0; i < N; ++i) {
      index = (index << 1) | bits[i];
    }
    return table[index];
  };
};
} // namespace functor
