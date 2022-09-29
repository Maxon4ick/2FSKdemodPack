cmake_minimum_required(VERSION 3.5)

project(anotherFun LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(anotherFun anotherFun.h)
target_link_libraries(${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
include("unittestA.cmake")
install(TARGETS anotherFun
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
