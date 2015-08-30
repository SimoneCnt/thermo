#
# Build the Thermo executable and the small thermo library
#

# Set the required CMake version
cmake_minimum_required(VERSION 2.8.7)

# Set the project name
project (Thermo C)

# Fix for MacOSX
set(CMAKE_MACOSX_RPATH 1)

# Compiling flags
set(CMAKE_C_FLAGS "-g -O3 -march=native -flto")

# Set the current directory for -I
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Source files
set(SOURCE_THERMO 
    src/calctermo.c 
    src/cumulvib.c 
    src/delete.c 
    src/diffthermo.c 
    src/init.c 
    src/printconfig.c 
    src/printthermo.c 
    src/readthermo.c 
    src/vdos.c 
)

# Build static libthermo.a
add_library(thermo-static STATIC ${SOURCE_THERMO})
set_target_properties(thermo-static PROPERTIES OUTPUT_NAME "thermo")
target_link_libraries(thermo-static m )

# Build shared libthermo.so
add_library(thermo-shared SHARED ${SOURCE_THERMO})
set_target_properties(thermo-shared PROPERTIES OUTPUT_NAME "thermo")
target_link_libraries(thermo-shared m )

# Build thermo (statically linked to libthermo)
add_executable(thermo src/thermo.c)
target_link_libraries(thermo thermo-static) 


