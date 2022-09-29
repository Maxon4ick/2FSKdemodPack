cmake_minimum_required(VERSION 3.0)
project (filter)




add_library(filter filter.cpp filter.h)



target_link_libraries(${PROJECT_NAME})

include("unittestFil.cmake")

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
