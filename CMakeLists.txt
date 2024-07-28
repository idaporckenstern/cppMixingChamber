project(MixingChamber)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

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
			outlet.h
			WriteData.cpp
			WriteData.h)

add_executable(mixing_chamber ${SOURCES})