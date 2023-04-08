project(reklibpp)

option(KSEQPP_READ_BUILD_STATIC "Build static libs. Off by default" OFF)

include(FetchContent)
FetchContent_Declare(
  zlib
  QUIET
  SYSTEM
  GIT_REPOSITORY       https://github.com/madler/zlib
  GIT_TAG              v1.2.13
  GIT_SHALLOW          TRUE
)
add_library(ZLIB STATIC IMPORTED)
if(KSEQPP_READ_BUILD_STATIC)
  set_target_properties(ZLIB PROPERTIES IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/lib/libz.a")
else()
  set_target_properties(ZLIB PROPERTIES IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/lib/libz.so")
endif()
target_include_directories(ZLIB PUBLIC ${ZLIB_SOURCE_DIR} ${ZLIB_BINARY_DIR})
set(ZLIB_BUILD_EXAMPLES OFF)
FetchContent_MakeAvailable(zlib)

add_library(kseqpp_read INTERFACE)
target_include_directories(kseqpp_read INTERFACE ${PROJECT_SOURCE_DIR}/kseqpp_read)
target_link_libraries(kseqpp_read INTERFACE ZLIB)
add_dependencies(kseqpp_read zlib)

# Builds the testing program. We use googletest as a testing framework
option(KSEQPP_READ_BUILD_TESTS "Build the benchmarks" OFF)
if (KSEQPP_READ_BUILD_TESTS)
  include(FetchContent)
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
  target_link_libraries(test kseqpp_read test_lib)
  add_test(NAME test COMMAND test)
endif() # BUILD_TESTS

option(KSEQPP_READ_BUILD_BENCHMARKS "Build the benchmarks" OFF)
if (KSEQPP_READ_BUILD_BENCHMARKS)
  include(ExternalProject)
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
  target_link_libraries(benchmark kseqpp_read ${STATIC_FLAGS})
endif()
