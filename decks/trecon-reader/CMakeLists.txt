cmake_minimum_required(VERSION 3.0)
project (vpic-reader CXX)

find_package (deltafs CONFIG REQUIRED)

add_executable (reader reader.cxx)
target_link_libraries (reader deltafs)
