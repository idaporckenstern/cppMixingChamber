cmake_minimum_required(VERSION 3.22.1)

project(MixingChamber)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

set(SOURCES main.cpp
			#DataTransfer.cpp
			#DataTransfer.h
			DataVector.cpp
			DataVector.h
			Functions.cpp
			Functions.h
			Mesh.cpp
			Mesh.h
			Opening.cpp
			Opening.h
			Outlet.cpp
			Outlet.h
			WriteData.cpp
			WriteData.h)

add_executable(mixing_chamber ${SOURCES})