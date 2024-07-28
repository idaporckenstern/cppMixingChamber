cmake_minimum_required(VERSION 3.22.1)

project(MixingChamber)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

add_library(data_vector STATIC DataVector.cpp)
add_library(functions STATIC functions.cpp)
add_library(mesh STATIC Mesh.cpp)
add_library(opening STATIC Opening.cpp)
add_library(outlet STATIC Outlet.cpp)
add_library(write_data STATIC WriteData.cpp)

# set(SOURCES main.cpp
# 			#DataTransfer.cpp
# 			#DataTransfer.h
# 			DataVector.cpp
# 			DataVector.h
# 			Functions.cpp
# 			Functions.h
# 			Mesh.cpp
# 			Mesh.h
# 			Opening.cpp
# 			Opening.h
# 			Outlet.cpp
# 			Outlet.h
# 			WriteData.cpp
# 			WriteData.h)

add_executable(mixing_chamber main.cpp)

target_link_libraries(mixing_chamber data_vector)
target_link_libraries(mixing_chamber functions)
target_link_libraries(mixing_chamber mesh)
target_link_libraries(mixing_chamber opening)
target_link_libraries(mixing_chamber outlet)
target_link_libraries(mixing_chamber write_data)