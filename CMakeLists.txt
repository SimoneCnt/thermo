#
# Thermo
#

cmake_minimum_required(VERSION 3.7.1)
set(CMAKE_MACOSX_RPATH 1)
project(Thermo C)
include(CheckIncludeFiles)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Set compiler options for C 
set(CMAKE_C_STANDARD 99)
# GCC
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU") 
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O1 -pedantic -Wall -W -Wextra -Wfatal-errors -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs ")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fshort-enums -fno-common -O3 -march=native -ffast-math")
# Clang
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O1 -pedantic -Wall -W -Wextra -Wfatal-errors -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fshort-enums -fno-common -O3 -march=native -ffast-math")
endif()

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
target_link_libraries(thermo m)

# Enable Testing
enable_testing()
add_subdirectory(examples EXCLUDE_FROM_ALL)

# Add debug and release targets
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

