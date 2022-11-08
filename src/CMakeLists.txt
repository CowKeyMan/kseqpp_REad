include(ExternalProject)

project(test)

# Builds the testing program. We use googletest as a testing framework
set(KSEQPP_READ_BUILD_TESTS OFF BOOL)
if (KSEQPP_READ_BUILD_TESTS)

find_package(ZLIB)

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
target_link_libraries(test ZLIB::ZLIB test_lib)
add_test(NAME test COMMAND test)

endif() # BUILD_TESTS
