//#include "OperatorsForMesons.h"

#include <benchmark/benchmark.h>

#include <cmath>

static void BM_exp(benchmark::State &state) {
  for (auto _ : state) {
    auto volatile x = std::exp(1.0);
  }
}

BENCHMARK(BM_exp);

#if 0
static void BM_create_momenta_old(benchmark::State &state) {
  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup{
    {0, {1, 2, 3}, {std::make_pair('>', 'x')}}
  };
}

BENCHMARK(BM_create_momenta_old);
#endif

BENCHMARK_MAIN();
