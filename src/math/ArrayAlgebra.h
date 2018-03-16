#pragma once

#include "interpolation/ValueTraits.h"
#include <array>
#include <type_traits>

template <class T, std::size_t N> struct ValueTraits<std::array<T, N>> {
  static constexpr std::array<T, N> zero() {
    std::array<T, N> out{};
    out.fill(0);
    return out;
  }
};

template <typename T, std::size_t N>
std::enable_if_t<std::is_arithmetic_v<T>, std::array<T, N> &>
operator+=(std::array<T, N> &a, std::array<T, N> const &b) {
  for (std::size_t i = 0; i < N; i++) {
    a[i] += b[i];
  }
  return a;
}

template <typename T, std::size_t N>
std::enable_if_t<std::is_arithmetic_v<T>, std::array<T, N> &>
operator-=(std::array<T, N> &a, std::array<T, N> const &b) {
  for (std::size_t i = 0; i < N; i++) {
    a[i] -= b[i];
  }
  return a;
}

template <typename T, std::size_t N>
std::enable_if_t<std::is_arithmetic_v<T>, std::array<T, N>>
operator+(std::array<T, N> const &a, std::array<T, N> const &b) {
  auto out = a;
  out += b;
  return out;
}

template <typename T, std::size_t N>
std::enable_if_t<std::is_arithmetic_v<T>, std::array<T, N>>
operator-(std::array<T, N> const &a, std::array<T, N> const &b) {
  auto out = a;
  out -= b;
  return out;
}

template <typename T, std::size_t N>
std::enable_if_t<std::is_arithmetic_v<T>, T>
operator*(std::array<T, N> const &a, std::array<T, N> const &b) {
  T dot = 0;
  for (std::size_t i = 0; i < N; i++) {
    dot += a[i] * b[i];
  }
  return dot;
}

template <typename T, typename S, std::size_t N>
std::enable_if_t<std::is_arithmetic_v<T> && std::is_assignable_v<T &, S>,
                 std::array<T, N> &>
operator*=(std::array<T, N> &a, S const &s) {
  for (std::size_t i = 0; i < N; i++) {
    a[i] *= s;
  }
  return a;
}

template <typename T, typename S, std::size_t N>
std::enable_if_t<std::is_arithmetic_v<T> && std::is_assignable_v<T &, S>,
                 std::array<T, N>>
operator*(S s, std::array<T, N> const &a) {
  auto out = a;
  out *= s;
  return out;
}

template <typename T, typename S, std::size_t N>
std::enable_if_t<std::is_arithmetic_v<T> && std::is_assignable_v<T &, S>,
                 std::array<T, N>>
operator*(std::array<T, N> const &a, S s) {
  auto out = a;
  out *= s;
  return out;
}

template <typename T, std::size_t N>
std::enable_if_t<std::is_arithmetic_v<T>, decltype(sqrt(std::declval<T>()))>
norm(std::array<T, N> const &a) {
  return sqrt(a*a);
}

template <typename T, std::size_t N>
std::enable_if_t<std::is_arithmetic_v<T>, std::array<T, N>>
normalized(std::array<T, N> const &a) {
  return a * (1.0 / norm(a));
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<T>, std::array<T, 3>>
cross(std::array<T, 3> const &a, std::array<T, 3> const &b) {
  return {a[1] * b[2] - a[2] * b[1], -a[0] * b[2] + a[2] * b[0],
          a[0] * b[1] - a[1] * b[0]};
}
