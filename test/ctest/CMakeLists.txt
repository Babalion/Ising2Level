cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(BUILD_ARCH "-m64")


add_executable(ctest_test_ising testIsing.cpp ../../SpinLattice2level.cpp)
add_test(ctest_test_ising ctest_exe_testIsing)
set_target_properties(ctest_test_ising PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")