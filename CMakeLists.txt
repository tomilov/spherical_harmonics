project(spherical_harmonics)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE ON)

cmake_minimum_required(VERSION 3.4 FATAL_ERROR)

include_directories(".")
include_directories("Math/")

add_executable(${PROJECT_NAME} "main.cpp" "Platform.cpp" "Math/Vector.cpp")

add_compile_options(-W -Wall)
set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS YES
)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++1z")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")

set(CMAKE_CXX_FLAGS_DEBUG "-ggdb -fno-inline -fno-omit-frame-pointer -DDEBUG=1 -D_GLIBCXX_DEBUG=1 -march=x86-64 -mtune=generic")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG=1 -march=native")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -gline-tables-only -fno-omit-frame-pointer -DNDEBUG=1 -march=native")
elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG=1 -march=native")
else()
    message(STATUS "Wrong build type selected, default to Debug.")
    set(CMAKE_BUILD_TYPE "Debug")
endif()
