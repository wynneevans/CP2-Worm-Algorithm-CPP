# set minium required version of CMake
cmake_minimum_required(VERSION 3.2)

# set the project name
project(CPP_code)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")
find_package(Boost 1.58.0 REQUIRED COMPONENTS system)

# include Boost library
include_directories(${Boost_INCLUDE_DIR})

# add the executable
add_executable(wormCP2Mu2dAFMplaq wormCP2Mu2dAFMplaq.cpp)

# link Boost library
target_link_libraries(wormCP2Mu2dAFMplaq LINK_PUBLIC ${Boost_LIBRARIES})
