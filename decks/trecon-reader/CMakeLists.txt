cmake_minimum_required(VERSION 3.0)
project (vpic-reader CXX)

find_package (MPI MODULE REQUIRED)
string (REPLACE " " ";" MPI_CXX_COMPILE_FLAGS_LIST "${MPI_CXX_COMPILE_FLAGS}")

#find_package (deltafs CONFIG REQUIRED)

# link shared lib with full rpath
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable (vpic-reader vpic-reader.cxx)

foreach (lcv ${MPI_CXX_COMPILE_FLAGS_LIST})
    if (NOT ${lcv} STREQUAL "")
        target_compile_options (vpic-reader PUBLIC ${lcv})
    endif ()
endforeach ()

# XXX: have to do this one dir at a time otherwise, otherwise I get
# error: target 'deltafs' INTERFACE_INCLUDE_DIRECTORIES contains path
#               prefixed in the source directory
foreach (lcv ${MPI_CXX_INCLUDE_PATH})
    target_include_directories (vpic-reader PUBLIC ${lcv})
endforeach ()

foreach (lcv ${MPI_CXX_LIBRARIES})
    target_link_libraries (vpic-reader ${lcv})
endforeach ()

set_property (TARGET vpic-reader APPEND PROPERTY LINK_FLAGS
              ${MPI_CXX_LINK_FLAGS})
#
# "make install" rules
#
install (TARGETS vpic-reader
         RUNTIME DESTINATION bin)

