find_package(GTest)

if(GTEST_FOUND)
  include_directories(${GTEST_INCLUDE_DIRS})

  add_executable(main
    main.cpp
    )

  target_link_libraries(main ${GTEST_LIBRARIES})

  enable_testing()
  add_test(NAME gtest
    COMMAND ./main)
else()
  message(STATUS "gtest was not found, building without this test.")
endif()
