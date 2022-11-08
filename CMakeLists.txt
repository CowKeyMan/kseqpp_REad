# This is the root cmake file where the sources lie
# Different CMake versions have more features.
# Since we will be using CUDA, we want a pretty high version

cmake_minimum_required(VERSION 3.19)

project(
  "SBWT Search"
  VERSION 0.0.1
  DESCRIPTION "A rewrite of kseqpp to only read sequences, and only read until a certain amount of characters, making it more scalable"
  HOMEPAGE_URL https://github.com/CowKeyMan/kseqpp_REad
  LANGUAGES CXX
)
set (CMAKE_CXX_STANDARD 17)

include(CheckLanguage)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(cmake/PreventBuildingInCmakeDirectory.cmake)
include(cmake/ExportCompileCommands.cmake)
include(cmake/CompilerWarnings.cmake)
include(cmake/Profiling.cmake)
include(cmake/Options.cmake)

add_subdirectory(src)
