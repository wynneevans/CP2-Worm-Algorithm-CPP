cmake_minimum_required(VERSION 3.19.1)
project(CP2_CPP_CODE CXX)

set(CMAKE_PREFIX_PATH "/Users/wynne/tools")
#find_package(Boost 1.74.0 REQUIRED COMPONENTS
#             Boost_format Boost_multi_array)

find_package(Boost 1.74.0 REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  add_executable(wormCP2Mu2dAFMplaq wormCP2Mu2dAFMplaq.cpp)
endif()

#target_link_libraries(wormCP2Mu2dAFMplaq Boost::format Boost::multi_array)


