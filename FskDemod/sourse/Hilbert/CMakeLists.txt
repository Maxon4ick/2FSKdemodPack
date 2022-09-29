cmake_minimum_required(VERSION 3.5)

project(Hilbert LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/include)
link_directories(${PROJECT_SOURCE_DIR}/lib)

add_library(Hilbert hilbert.cpp hilbert.h)

target_link_libraries(Hilbert libfftw3.a)
target_link_libraries(${PROJECT_NAME} filter)
include("unittestHil.cmake")
install(TARGETS Hilbert
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
