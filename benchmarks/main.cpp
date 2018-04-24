#include "OperatorsForMesons.h"

#include <benchmark/benchmark.h>

#include <cmath>

static void BM_exp(benchmark::State &state) {
  for (auto _ : state) {
    auto volatile x = std::exp(1.0);
  }
}

BENCHMARK(BM_exp);

static void BM_create_momenta_old(benchmark::State &state) {
  ssize_t Lx = 32, Ly = 32, Lz = 32, Lt = 96;
  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup{{0, {1, 2, 3}, {}}};
  array_cd_d2 momentum(boost::extents[1][Lx*Ly*Lz]);

  for (auto _ : state) {
    create_momenta(Lx, Ly, Lz, vdaggerv_lookup, momentum);
  }
}

BENCHMARK(BM_create_momenta_old);

BENCHMARK_MAIN();
