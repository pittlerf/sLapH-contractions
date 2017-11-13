// Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

#pragma once

#include <cassert>
#include <iomanip>
#include <map>
#include <iostream>
#include <stdexcept>
#include <type_traits>

enum class DilutionType { block, interlace };

std::map<DilutionType const, std::string const> dilution_names = {
    {DilutionType::block, "B"}, {DilutionType::interlace, "I"}};

class BlockIterator {
public:
  struct Element {
    int source;
    int sink;
  };

  BlockIterator(int const slice_source,
                int const slice_sink,
                int const block_source,
                int const block_sink,
                int const num_slice,
                int const num_block,
                int const pass,
                DilutionType const type)
      : slice_source_(slice_source),
        slice_sink_(slice_sink),
        block_source_(block_source),
        block_sink_(block_sink),
        num_slice_(num_slice),
        num_block_(num_block),
        pass_(pass),
        type_(type) {}

  Element operator*() { return {slice_source_, slice_sink_}; }

  BlockIterator operator++() {
    if (type_ == DilutionType::block) {
      int const block_size = num_slice_ / num_block_;
      int const block_source_begin = block_source_ * block_size;
      int const block_source_end = (block_source_ + 1) * block_size;
      int const block_sink_begin = block_sink_ * block_size;
      int const block_sink_end = (block_sink_ + 1) * block_size;

      // Go to the next sink slice.
      ++slice_sink_;

      // If we iterated through all sinks in this block, we need to go to the
      // next source.
      if (slice_sink_ == block_sink_end) {
        slice_sink_ = block_sink_begin;
        ++slice_source_;
      }

      // If we iterated through all the sources in our block, we need to
      // exchange
      // source and sink block and repeat.
      if (slice_source_ == block_source_end) {
        slice_source_ = block_sink_begin;
        slice_sink_ = block_source_begin;
        ++pass_;

        if (block_source_ == block_sink_) {
          ++pass_;
        }

        std::swap(block_source_, block_sink_);
      }
    } else if (type_ == DilutionType::interlace) {
      int const block_size = num_slice_ / num_block_;
      int const block_source_begin = block_source_;
      int const block_sink_begin = block_sink_;

      // Go to the next sink slice.
      slice_sink_ += block_size;

      // If we iterated through all sinks in this block, we need to go to the
      // next source.
      if (slice_sink_ >= num_slice_) {
        slice_sink_ = block_sink_begin;
        slice_source_ += block_size;
      }

      // If we iterated through all the sources in our block, we need to
      // exchange
      // source and sink block and repeat.
      if (slice_source_ >= num_slice_) {
        slice_source_ = block_sink_begin;
        slice_sink_ = block_source_begin;
        ++pass_;

        if (block_source_ == block_sink_) {
          ++pass_;
        }

        std::swap(block_source_, block_sink_);
      }
    } else {
      throw std::domain_error("This dilution scheme is not implemented.");
    }

    return *this;
  }

  bool operator!=(BlockIterator const &other) const {
    return block_source_ != other.block_source_ ||
           block_sink_ != other.block_sink_ ||
           num_slice_ != other.num_slice_ || num_block_ != other.num_block_ ||
           pass_ != other.pass_ || type_ != other.type_;
  }

private:
  int slice_source_;
  int slice_sink_;
  int block_source_;
  int block_sink_;
  int num_slice_;
  int num_block_;
  int pass_;
  DilutionType type_;
};

class DilutionIterator {
public:
  struct Element {
    int source;
    int sink;
  };

  DilutionIterator(int const block_source,
                   int const block_sink,
                   int const num_slice,
                   int const num_block,
                   DilutionType const type)
      : block_source_(block_source),
        block_sink_(block_sink),
        num_slice_(num_slice),
        num_block_(num_block),
        type_(type) {}

  DilutionIterator operator*() const { return *this; }

  DilutionIterator operator++() {
    ++block_sink_;

    if (block_sink_ == num_block_) {
      ++block_source_;
      block_sink_ = block_source_;
    }

    return *this;
  }

  bool operator!=(DilutionIterator const &other) const {
    return block_source_ != other.block_source_ ||
           block_sink_ != other.block_sink_ ||
           num_slice_ != other.num_slice_ || num_block_ != other.num_block_;
  }

  BlockIterator begin() const {
    int const block_size = num_slice_ / num_block_;

    if (type_ == DilutionType::block) {
      return BlockIterator(block_source_ * block_size,
                           block_sink_ * block_size,
                           block_source_,
                           block_sink_,
                           num_slice_,
                           num_block_,
                           0,
                           type_);
    } else {
      return BlockIterator(block_source_,
                           block_sink_,
                           block_source_,
                           block_sink_,
                           num_slice_,
                           num_block_,
                           0,
                           type_);
    }
  }

  BlockIterator end() const {
    int const block_size = num_slice_ / num_block_;

    if (type_ == DilutionType::block) {
      return BlockIterator(block_source_ * block_size,
                           block_sink_ * block_size,
                           block_source_,
                           block_sink_,
                           num_slice_,
                           num_block_,
                           2,
                           type_);
    } else {
      return BlockIterator(block_source_,
                           block_sink_,
                           block_source_,
                           block_sink_,
                           num_slice_,
                           num_block_,
                           2,
                           type_);
    }
  }

  int source() const { return block_source_; }

  int sink() const { return block_sink_; }

private:
  int block_source_;
  int block_sink_;
  int num_slice_;
  int num_block_;
  DilutionType type_;
};

class DilutionScheme {
public:
  DilutionScheme(int const num_slice,
                 int const num_block,
                 DilutionType const type)
      : num_slice_(num_slice), num_block_(num_block), type_(type) {}

  DilutionIterator operator[](int const i) const {
    int block_sink = 0;
    int block_source = 0;
    for (int j = 0; j < i; ++j) {
      ++block_sink;

      if (block_sink == num_block_) {
        ++block_source;
        block_sink = block_source;
      }
    }

    return DilutionIterator(
        block_source, block_sink, num_slice_, num_block_, type_);
  }

  int size() const { return num_block_ * (num_block_ + 1) / 2; }

  DilutionIterator begin() const {
    return DilutionIterator(0, 0, num_slice_, num_block_, type_);
  }

  DilutionIterator end() const {
    return DilutionIterator(
        num_block_, num_block_, num_slice_, num_block_, type_);
  }

 private:
  int num_slice_;
  int num_block_;
  DilutionType type_;
};

inline void test_dilution_scheme(int const num_slice, int const num_block, DilutionType const type) {
  auto const name = dilution_names[type];
  std::cout << "T = " << num_slice << ", T" << name
            << num_block << " (Morningstar), T" << name << (num_slice / num_block)
            << " (Other):\n\n";

  DilutionScheme dilution_scheme(num_slice, num_block, type);
  for (int b = 0; b < dilution_scheme.size(); ++b) {
    auto const blocks = dilution_scheme[b];
    std::cout << std::setw(2) << blocks.source() << " => " << std::setw(2)
              << blocks.sink() << "\n";

    for (auto const slices : blocks) {
      std::cout << "  " << std::setw(2) << slices.source << " -> "
                << std::setw(2) << slices.sink << "\n";
    }
  }

  std::cout << "\n\n";
}
