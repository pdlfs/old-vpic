cmake_minimum_required(VERSION 3.0)
project (vpic-reader CXX)

find_package (deltafs CONFIG REQUIRED)

# link shared lib with full rpath
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable (vpic-reader vpic-reader.cxx)
target_link_libraries (vpic-reader deltafs)

#
# "make install" rules
#
install (TARGETS vpic-reader
         RUNTIME DESTINATION bin)

