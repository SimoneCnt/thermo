#
# Build the Thermo executable and the small thermo library
#

# Set the required CMake version
cmake_minimum_required(VERSION 2.8.7)

# Set the project name
project (Thermo C)

# Fix for MacOSX
set(CMAKE_MACOSX_RPATH 1)

# Set the current directory for -I
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Source files
set(SOURCE_THERMO calctermo.c cumulvib.c delete.c diffthermo.c init.c printconfig.c printthermo.c readthermo.c vdos.c vdosfvib.c)

# Build static libthermo.a
add_library(thermo-static STATIC ${SOURCE_THERMO})
set_target_properties(thermo-static PROPERTIES OUTPUT_NAME "thermo")
target_link_libraries(thermo-static m )

# Build shared libthermo.so
add_library(thermo-shared SHARED ${SOURCE_THERMO})
set_target_properties(thermo-shared PROPERTIES OUTPUT_NAME "thermo")
target_link_libraries(thermo-shared m )


# Build thermo (statically linked to libthermo)
add_executable(thermo thermo.c)
target_link_libraries(thermo thermo-static) 


