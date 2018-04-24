#include <benchmark/benchmark.h>

#include <cmath>

static void BM_exp(benchmark::State &state) {
  for (auto _ : state) {
    auto volatile x = std::exp(1.0);
  }
}

BENCHMARK(BM_exp);

BENCHMARK_MAIN();
