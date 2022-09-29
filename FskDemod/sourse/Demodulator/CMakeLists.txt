cmake_minimum_required(VERSION 3.5)

project(Demodulator LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE  SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}*.cpp SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}*.h)


add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} Hilbert)
target_link_libraries(${PROJECT_NAME} filter)

