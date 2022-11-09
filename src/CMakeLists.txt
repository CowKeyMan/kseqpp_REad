project(kseqpp_read)
include_directories(${PROJECT_SOURCE_DIR}/kseqpp_read)
message(${PROJECT_SOURCE_DIR}/kseqpp_read)

# Builds the testing program. We use googletest as a testing framework
set(KSEQPP_READ_BUILD_TESTS OFF BOOL)
if (KSEQPP_READ_BUILD_TESTS)
  include(FetchContent)
  find_package(ZLIB)
  FetchContent_Declare(
    googletest
    QUIET
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        release-1.12.1
  )
  FetchContent_MakeAvailable(googletest)
  enable_testing()
  add_library(test_lib INTERFACE)
  target_link_libraries(
    test_lib INTERFACE
    gtest
    gtest_main
    gmock
    gmock_main
    gcov
  )

  add_executable(
    test
    "${PROJECT_SOURCE_DIR}/test.cpp"
  )
  target_link_libraries(test ZLIB::ZLIB test_lib)
  add_test(NAME test COMMAND test)
endif() # BUILD_TESTS

set(KSEQPP_READ_BUILD_BENCHMARKS OFF BOOL)
if (KSEQPP_READ_BUILD_BENCHMARKS)
  include(ExternalProject)
  find_package(ZLIB)
  if(EXISTS "${CMAKE_BINARY_DIR}/external/kseqpp/lib/pkgconfig/kseq++.pc")
    set(KSEQPP_FOUND TRUE)
  else()
    set(KSEQPP_FOUND FALSE)
  endif()
  if (NOT KSEQPP_FOUND)
    ## Fetch kseqpp
    ExternalProject_Add(
      kseqpp
      GIT_REPOSITORY https://github.com/cartoonist/kseqpp
      GIT_TAG        v0.2.1
      PREFIX         "${CMAKE_BINARY_DIR}/external/kseqpp"
      CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    )
  endif()
  include_directories(SYSTEM "${CMAKE_BINARY_DIR}/external/kseqpp/include")
  add_executable(
    benchmark
    "${PROJECT_SOURCE_DIR}/benchmark.cpp"
  )
  if(NOT KSEQPP_FOUND)
    add_dependencies(benchmark kseqpp)
  endif()
  target_link_libraries(benchmark ZLIB::ZLIB)
endif()
