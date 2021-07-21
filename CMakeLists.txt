cmake_minimum_required(VERSION 3.19)
project(Ising2Level)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pthread")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -g -O0")
set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -g -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=[sanitizer_name] [additional_options] [-g] [-OX]")

set(CMAKE_C_FLAGS_RELEASE "-Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set(BUILD_ARCH "-m64")

#######################################################################################################
#comment out the following lines if you dont need plots
find_package(OpenCV REQUIRED)
add_subdirectory(lib/cv-plot-1.2.1/CvPlot)

add_executable(ising-with-plots ising-with-plots.cpp SpinLattice2level.cpp)
target_link_libraries(ising-with-plots ${OpenCV_LIBS} CvPlot stdc++fs)
set_target_properties(ising-with-plots PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(ising-live ising-live.cpp SpinLattice2level.cpp)
target_link_libraries(ising-live ${OpenCV_LIBS} CvPlot stdc++fs)
set_target_properties(ising-live PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

#######################################################################################################
#with included plotting tools
add_executable(ising-headless ising-headless.cpp SpinLattice2level.cpp)
target_link_libraries(ising-headless stdc++fs)
set_target_properties(ising-headless PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

enable_testing()
add_subdirectory(test/ctest)