cmake_minimum_required(VERSION 3.0)
project (vpic-reader CXX)

find_package (deltafs CONFIG REQUIRED)

# link shared lib with full rpath
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable (reader reader.cxx)
target_link_libraries (reader deltafs)

#
# "make install" rules
#
install (TARGETS reader
         RUNTIME DESTINATION bin)

