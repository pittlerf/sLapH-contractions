#include "ComplexProduct.hpp"

ComplexProduct operator+(ComplexProduct left, ComplexProduct const &right) {
  left += right;
  return left;
}

ComplexProduct operator*(ComplexProduct const &left, ComplexProduct const &right) {
  return {left.full * right.full, left.select * right.select};
}

ComplexProduct operator/(ComplexProduct const &left, double const &right) {
  return {left.full / right, left.select / right};
}

ComplexProduct make_complex_product(Complex const &other, bool const ignore_imag) {
  ComplexProduct rval{other, other};
  if (ignore_imag) {
    rval.select.imag(0.0);
  }
  return rval;
}
