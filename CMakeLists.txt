cmake_minimum_required(VERSION 3.16)
project(mptscfl CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(MPTSCFL_WITH_GUROBI "Build Gurobi-dependent targets" ON)

# ---- core (header-only) + tests: no solver required -------------------------
add_executable(test_core tests/test_core.cpp)
enable_testing()
# Pass an instance path to test_core to also exercise the parser on real data.
add_test(NAME core COMMAND test_core)

# ---- Gurobi ------------------------------------------------------------------
# Set GUROBI_HOME (e.g. C:/gurobi1300/win64 or /opt/gurobi1300/linux64) or pass
# -DGUROBI_HOME=... ; links gurobi_c++ + the versioned C library.
if(MPTSCFL_WITH_GUROBI)
  if(NOT DEFINED GUROBI_HOME)
    set(GUROBI_HOME $ENV{GUROBI_HOME})
  endif()
  find_path(GUROBI_INCLUDE_DIR gurobi_c++.h HINTS ${GUROBI_HOME}/include)
  file(GLOB GUROBI_C_LIB_CANDIDATES
       ${GUROBI_HOME}/lib/libgurobi*.so ${GUROBI_HOME}/lib/gurobi*.lib)
  find_library(GUROBI_CXX_LIB NAMES gurobi_c++ gurobi_c++md2017
               HINTS ${GUROBI_HOME}/lib)
  if(GUROBI_INCLUDE_DIR AND GUROBI_CXX_LIB AND GUROBI_C_LIB_CANDIDATES)
    list(GET GUROBI_C_LIB_CANDIDATES 0 GUROBI_C_LIB)
    message(STATUS "Gurobi: ${GUROBI_INCLUDE_DIR} | ${GUROBI_CXX_LIB} | ${GUROBI_C_LIB}")
    add_executable(mptscfl
      src/main.cpp src/exact_model.cpp src/single_stage_model.cpp
      src/two_steps_solver.cpp src/benders_model.cpp src/lagrangian.cpp)
    target_include_directories(mptscfl PRIVATE ${GUROBI_INCLUDE_DIR})
    target_link_libraries(mptscfl PRIVATE ${GUROBI_CXX_LIB} ${GUROBI_C_LIB})
  else()
    message(WARNING "Gurobi not found (set GUROBI_HOME); building only test_core.")
  endif()
endif()
