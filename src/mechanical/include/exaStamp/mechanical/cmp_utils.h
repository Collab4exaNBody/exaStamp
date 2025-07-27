#pragma once

#include <cmath>

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

inline constexpr size_t upper_triangle_index(size_t i, size_t j) {
  size_t ii = cmp::min(i, j);
  size_t jj = cmp::max(i, j);
  return jj * (jj + 1) / 2 + ii;
}

inline constexpr void unflatten_lower(size_t flat_index, size_t& i, size_t& j) {
  i = static_cast<size_t>(cmp::floor((std::sqrt(1 + 8 * flat_index) - 1) / 2));
  j = flat_index - i * (i + 1) / 2;
}

} // namespace cmp
