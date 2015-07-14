cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS "-stdlib=libstdc++")
project( fhis )
find_package(OpenCV REQUIRED)
find_package(CGAL QUIET COMPONENTS Core )

#include(src)
include( ${CGAL_USE_FILE} )
include( CGAL_CreateSingleSourceCGALProgram )

set(EXECUTABLE_OUTPUT_PATH ../bin)

include_directories( src )
include_directories( ${OpenCV_INCLUDE_DIRS} )

file(GLOB_RECURSE fhis_SRC "src/*.cpp")
add_executable( fhis ${fhis_SRC})
target_link_libraries( fhis ${OpenCV_LIBS})
