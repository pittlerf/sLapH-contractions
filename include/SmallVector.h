#pragma once

#include <array>
#include <cassert>

/*! Vector with fixed capacity.

  The `std::vector` uses a dynamic allocation, this might cost us performance when we do
  not need an arbitrary amount of data. In some parts of the code, the needed capacity
  will only be around 10. Therefore it is sufficient to use `std::array` which is
  completely allocated on the stack and therefore without `malloc` or `new` and also
  without indirections.

  This class uses `std::array` internally but provides the distinction between “size” and
  “capacity” that `std::vector` gives us.

  In contrast to the `SmallVector` of the LLVM project this one does not fall back to a
  dynamic allocation if the capacity is reached. Therefore it is undefined to write beyond
  the capacity. The @ref push_back member function contains an `assert` in debugging mode.

  @tparam Data Type of the payload
  @tparam capacity_ Fixed capacity and maximal size of the container.
  */
template <typename Data, int capacity_>
class SmallVector {
  public:
    SmallVector() {}

    SmallVector(std::vector<Data> const &other) {
      if (other.size() > capacity_) {
        throw std::runtime_error("Size of other vector exceeds SmallVector capacity");
      }

      std::copy(std::begin(other), std::end(other), std::begin(data_));
      size_ = other.size();
    }

    int size() const { return size_; }

    int capacity() const { return capacity_; }

    Data *begin() const { return &data_[0]; }

    Data *end() const { return &data_[size_ + 1]; }

    Data &operator[](int const i) { return data_[i]; }

    Data &at(int const i) {
      if (i < 0 || size_ < i) {
        throw std::runtime_error("Index outside of size.");
      }

      return data_[i];
    }

    void push_back(Data const &elem) {
      assert(size_ < capacity_ && "SmallVector must not be full when using push_back.");
      data_[size_] = elem;
      ++size_;
    }

   private:
    std::array<Data, capacity> data_;

    int size_ = 0;
}
