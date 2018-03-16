#pragma once

#include <eigen3/Eigen/Dense>

template <class T> struct ValueTraits {

  static constexpr T zero() {
    if constexpr (std::is_base_of_v<Eigen::DenseBase<T>, T>) {
      return T::Zero();
    } else {
      return 0;
    }
  }
};
