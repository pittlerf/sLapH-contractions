find_package(benchmark)

if(${benchmark_FOUND})
  add_executable(bench-main
    main.cpp

    accumulation.cpp
    create_momenta.cpp
    omp_wtime.cpp

    )

  target_link_libraries(bench-main benchmark::benchmark lcontract)

else()
  message(STATUS "Google benchmark framework was not found, not building benchmarks.")
endif()
