#pragma once

#include <tuple>
#include <utility>

template <class U, class... T> constexpr bool is_all_same() {
  return (std::is_same_v<U, T> && ... && true);
}

template <int I, typename... Ts>
using getType = typename std::tuple_element<I, std::tuple<Ts...>>::type;

template <int I, class... Ts> constexpr decltype(auto) get(Ts &&... ts) {
  return std::get<I>(std::forward_as_tuple(ts...));
}

template <int First, int Last, class Lambda>
inline void static_for(Lambda const &f) {

  if constexpr (First < Last) {
    f(std::integral_constant<int, First>{});
    static_for<First + 1, Last>(f);
  }
}

template <typename> struct ApplyNTimes_impl;

template <std::size_t... I> struct ApplyNTimes_impl<std::index_sequence<I...>> {
  template <typename V, int N> using type = V;

  template <template <typename...> typename T, typename V>
  using apply = T<type<V, I>...>;
};

template <int N>
using ApplyNTimes = ApplyNTimes_impl<std::make_index_sequence<N>>;
