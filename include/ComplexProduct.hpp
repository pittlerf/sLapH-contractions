#pragma once

#include <complex>

typedef std::complex<double> Complex;

struct ComplexProduct {
  Complex full;
  Complex select;

  ComplexProduct &operator+=(ComplexProduct const &right) {
    full += right.full;
    select += right.select;

    return *this;
  }
};

ComplexProduct operator+(ComplexProduct const &left, ComplexProduct const &right);

ComplexProduct operator*(ComplexProduct const &left, ComplexProduct const &right);

ComplexProduct operator/(ComplexProduct const &left, double const &right);

ComplexProduct make_complex_product(Complex const &other, bool const ignore_imag);
