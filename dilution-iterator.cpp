// Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <type_traits>

struct Combination {
  int block;
  int source;
  int sink;
};

class DilutionIterator {
public:
  DilutionIterator(int const block,
                   int const source,
                   int const sink,
                   int const num_slice,
                   int const num_block)
      : block_(block),
        source_(source),
        sink_(sink),
        num_slice_(num_slice),
        num_block_(num_block) {}

  Combination operator*() const { return {block_, source_, sink_}; }

  bool operator!=(DilutionIterator const &other) {
    return block_ != other.block_ || source_ != other.source_ ||
           sink_ != other.sink_ || num_slice_ != other.num_slice_ ||
           num_block_ != other.num_block_;
  }

protected:
  int block_;
  int source_;
  int sink_;
  int num_slice_;
  int num_block_;
};

class BlockDilutionIterator : public DilutionIterator {
public:
  using DilutionIterator::DilutionIterator;

  BlockDilutionIterator operator++() {
    int const block_size = num_slice_ / num_block_;
    int const block_begin = block_ * block_size;
    int const block_end = (block_ + 1) * block_size;

    // Go to the next sink.
    ++sink_;

    // If we iterated through all sinks, we must go to the next source in
    // our block.
    if (sink_ == block_end) {
      sink_ = block_begin;
      ++source_;
    }

    // If we iterated through all the sources, we must go to the next
    // block.
    if (source_ == block_end) {
      source_ = block_end;
      sink_ = block_end;
      ++block_;
    }

    return *this;
  }
};

class InterlaceDilutionIterator : public DilutionIterator {
public:
  using DilutionIterator::DilutionIterator;

  InterlaceDilutionIterator operator++() {
    // Go to the next sink.
    sink_ += num_block_;

    // If we iterated through all sinks, we must go to the next source in
    // our block.
    if (sink_ >= num_slice_) {
      sink_ = block_;
      source_ += num_block_;
    }

    // If we iterated through all the sources, we must go to the next
    // block.
    if (source_ >= num_slice_) {
      ++block_;
      sink_ = block_;
      source_ = block_;
    }

    return *this;
  }
};

template <typename Iterator>
class Iteration {
public:
  Iteration(int const num_slice, int const num_block)
      : num_slice_(num_slice), num_block_(num_block) {}

  Iterator begin() const { return Iterator{0, 0, 0, num_slice_, num_block_}; }

  Iterator end() const {
    if (std::is_same<Iterator, BlockDilutionIterator>::value) {
      return Iterator{num_block_, num_slice_, num_slice_, num_slice_,
                      num_block_};
    }
    else if (std::is_same<Iterator, InterlaceDilutionIterator>::value) {
      return Iterator{num_block_, num_block_, num_block_, num_slice_,
                      num_block_};
    }
  }

private:
  int num_slice_;
  int num_block_;
};

int main() {
  std::cout << "Block:\n\n";
  for (auto comb : Iteration<BlockDilutionIterator>(48, 24)) {
    std::cout << std::setw(2) << comb.block << ": " << std::setw(2)
              << comb.source << " -> " << std::setw(2) << comb.sink << "\n";
  }

  std::cout << "\n\nInterlace:\n\n";

  for (auto comb : Iteration<InterlaceDilutionIterator>(8, 2)) {
    std::cout << std::setw(2) << comb.block << ": " << std::setw(2)
              << comb.source << " -> " << std::setw(2) << comb.sink << "\n";
  }
}
