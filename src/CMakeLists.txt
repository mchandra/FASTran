# To compile the code, make a seperate build folder and then do the following:
#
# cd /location/to/build
# cmake /location/to/grim/src
# make
#
# If a different compiler than the one cmake detects is needed, then the
# recommended way is to do the following:
#
# cmake -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpic++ /location/to/grim/src 

cmake_minimum_required(VERSION 2.8)

project(grim)


# ------------------------------USER OPTIONS----------------------------------#
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O3")

# Set custom install folders here
set(PETSC_DIR "/home/mc/data/mc/Downloads/petsc_optimized")

# Double or float precision?
set(REAL "double")

# Use as ARRAY_ARGS ARRAY_SIZE while declaring arrays in functions arguments
# Must use -std=c99 if using the "restrict" keyword while compiling 
set(ARRAY_ARGS "const restrict static")

set(DUMP_FILE_PREFIX "\"data\"")

# Time stepping options
set(DT "0.0001")
set(DT_DUMP "0.0001")
set(START_TIME "0.")
set(FINAL_TIME "4.")
set(COURANT "10000.")
# Choose between BACKWARD_EULER or CRANK_NICHOLSON
set(TIME_STEPPING SEMI_IMPLICIT)

# Domain size
set(N1    "200")
set(N2    "200")

set(X1_A  "-1.")
set(X1_B  "1.")
set(X2_A  "-1.")
set(X2_B  "1.")

# Limiter?
# Choose between MINMOD, VAN_LEER and MC
set(LIMITER "MC")

# Diffusion coefficient
set(D0 "1.")

# Initial condition parameters
set(INITIAL_RADIUS "0.05")
set(TEMP_INSIDE    "10")
set(TEMP_OUTSIDE   "0.01")

set(USE_OPENMP   "NO")


# ------------------------------END OF USER OPTIONS---------------------------#

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Need to find all packages before starting compilation, so that all headers
# needed by the code are included.

# Should the external libraries needed by the internal libraries go into their
# own CMakeLists in their own folders?

# TODO: Only find packages depending on what options are set.

# System libraries
find_library(MATH_LIBRARIES m REQUIRED)

# External packages
find_package(PETSc REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${PETSC_INCLUDES})

configure_file(config.h.in ${CMAKE_SOURCE_DIR}/inputs.h)

add_executable(implicitconduction implicitconduction.c implicitconduction.h
               timestepper.h timestepper.c diagnostics.c initialconditions.c
               residual.c computeflux.c limiters.h limiters.c)

target_link_libraries(implicitconduction
                      ${MATH_LIBRARIES}
                      ${PETSC_LIBRARIES})

message("")
message("#################")
message("# Build options #")
message("#################")
message("")
message("C Compiler       : " ${CMAKE_C_COMPILER})
message("C_FLAGS          : " ${CMAKE_C_FLAGS})
message("PETSc dir        : " ${PETSC_DIR})
message("")
