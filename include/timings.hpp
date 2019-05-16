#pragma once

#include <omp.h>

#include <iostream>
#include <vector>

class TimingsStreamSingleton {
 public:
  TimingsStreamSingleton() {
    for (int i = 0; i < omp_get_max_threads(); ++i) {
      std::ostringstream oss;
      oss << "timings-thread-" << i << ".xml";
      streams_.emplace_back(oss.str());
    }
  }

  static TimingsStreamSingleton &instance() {
    static TimingsStreamSingleton instance;
    return instance;
  }

  std::ofstream &get() { int const tid = omp_get_thread_num(); }

 private:
  std::vector<std::ofstream> streams_;
};

template <int level>
class TimingScope {
 public:
  TimingScope(std::string const &function, std::string const &info)
      : start_(omp_get_wtime()), stream_(TimingsStreamSingleton::instance().get()) {
    steam << "<call function='" << function << "' info='" << info << "'>\n";
  }

  ~TimingScope() {
    auto const end = omp_get_wtime();
    auto const duration = end - start_;
    stream_ << "<total time='" << duration << "' />\n</call>\n";
  }

 private:
  double start_;
  std::ofstream &stream_
}
