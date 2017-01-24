#
# Build the Thermo executable
#

# Set the required CMake version
cmake_minimum_required(VERSION 2.8.7)

# Set the project name
project (Thermo C)

# Fix for MacOSX
set(CMAKE_MACOSX_RPATH 1)

# Set basic compiler options
set(CMAKE_C_FLAGS "")
if((${CMAKE_C_COMPILER_ID} MATCHES "GNU") OR (${CMAKE_C_COMPILER_ID} MATCHES "Clang"))
    set(CMAKE_C_FLAGS "${MCAKE_C_FLAGS} -std=c99 -g -O2")
elseif(${CMAKE_C_COMPILER_ID} MATCHES "PGI")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99 -gopt")
endif()

# Set warnings option
if((${CMAKE_C_COMPILER_ID} MATCHES "GNU") OR (${CMAKE_C_COMPILER_ID} MATCHES "Clang"))
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pedantic -Wall -W -Wextra -Wfatal-errors -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs")
elseif(${CMAKE_C_COMPILER_ID} MATCHES "PGI")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ")
endif()

# Set optimizations
if((${CMAKE_C_COMPILER_ID} MATCHES "GNU") OR (${CMAKE_C_COMPILER_ID} MATCHES "Clang"))
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fshort-enums -fno-common -O3 -march=native -flto")
elseif(${CMAKE_C_COMPILER_ID} MATCHES "PGI")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ")
endif()

message("Compiling options: ${CMAKE_C_FLAGS}")

# Set the current directory for -I
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

# Source files
set(SRC_THERMO 
    src/calctermo.c 
    src/cumulvib.c 
    src/delete.c 
    src/diffthermo.c 
    src/init.c 
    src/printconfig.c 
    src/printthermo.c 
    src/readthermo.c 
    src/vdos.c 
    src/thermo.c
)

# Build thermo (statically linked to libthermo)
add_executable(thermo ${SRC_THERMO})

# Link math library
IF(UNIX)
  target_link_libraries(thermo m)
ENDIF(UNIX)

