cmake_minimum_required(VERSION 3.19)
project(Ising2Level)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pthread")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set(BUILD_ARCH "-m64")

#######################################################################################################
#comment out the following lines if you dont need plots
find_package(OpenCV REQUIRED)
add_subdirectory(lib/cv-plot-1.2.1/CvPlot)

add_executable(ising ising.cpp SpinLattice2level.cpp)
target_link_libraries(ising ${OpenCV_LIBS} CvPlot stdc++fs)
set_target_properties(ising PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

#######################################################################################################
#without included plotting tools
add_executable(isingHeadless isingHeadless.cpp SpinLattice2level.cpp)
target_link_libraries(isingHeadless stdc++fs)
set_target_properties(isingHeadless PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

enable_testing()
add_subdirectory(test/ctest)