cmake_minimum_required(VERSION 3.10)
project(toddcox-faster)

set(CMAKE_CXX_STANDARD 17)

option(TC_BUILD_EXAMPLE "Build example executables" OFF)

add_subdirectory(ext)

add_library(tc INTERFACE)
target_link_libraries(tc INTERFACE eigen)

target_include_directories(tc INTERFACE include)

if (TC_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()

if (TC_BUILD_TESTING)
    enable_testing()
    add_subdirectory(Testing)
endif()
