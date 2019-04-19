#
# Thermo
#

cmake_minimum_required(VERSION 3.7.1)
set(CMAKE_MACOSX_RPATH 1)
project(Thermo C)
include(CheckIncludeFiles)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Set compiler options for C 
set(CMAKE_C_STANDARD 11)
# GCC
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU") 
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O1 -pedantic -Wall -W -Wextra -Wfatal-errors -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs ")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fshort-enums -fno-common -O3 -march=native")
# Clang
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O1 -pedantic -Wall -W -Wextra -Wfatal-errors -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fshort-enums -fno-common -O3 -march=native")
endif()

# Get version from last git commit
execute_process(COMMAND git describe --long --tags --dirty --always  OUTPUT_VARIABLE GIT_VERSION RESULT_VARIABLE WITH_GIT OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT (${WITH_GIT} EQUAL 0))
    set(GIT_VERSION "unknown")
endif(NOT (${WITH_GIT} EQUAL 0))

# Set include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src)

# Generate thermo.h
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/thermo.h.in ${CMAKE_CURRENT_BINARY_DIR}/src/thermo.h)

# Source files
set(SRC_THERMO 
    src/getline.c
    src/calctermo.c 
    src/cumulvib.c 
    src/delete.c 
    src/diffthermo.c 
    src/init.c 
    src/printconfig.c 
    src/printthermo.c 
    src/readthermo.c 
    src/readhessian.c
    src/calcfreqs.c
    src/dsyev.c
    src/vdos.c 
    src/thermo.c
    src/kernel.c
    src/solvents.c
)

# Build thermo (statically linked to libthermo)
add_executable(thermo ${SRC_THERMO})

# Link math library
target_link_libraries(thermo m)

# Linear algebra (lapack and blas)
if(DEFINED ENV{LAPACK_INCL} AND DEFINED ENV{LAPACK_LINK})
    include_directories($ENV{LAPACK_INLC})
    target_link_libraries(thermo $ENV{LAPACK_LINK})
    add_definitions(-DHAVE_LAPACK)
else()
    find_package(LAPACK)
    find_package(BLAS)
    if(${LAPACK_FOUND} AND ${BLAS_FOUND})
        target_link_libraries(thermo ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        add_definitions(-DHAVE_LAPACK)
    else()
        message(WARNING "WARNING!\nImpossible to find lapack. Some functions will be disabled.\nSet environmental variables LAPACK_INCL and LAPACK_LINK")
    endif()
endif()

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

# Add a target to generate documentation with pandoc
add_custom_target(doc
    pandoc --filter pandoc-citeproc -s README.md INSTALL.md USAGE.md EXAMPLES.md THEORY.md doc/references.md doc/metadata.yaml -V documentclass=article -V classoption:fleqn -V geometry:margin=1.5in -V colorlinks=true --table-of-contents --number-sections -o ${CMAKE_SOURCE_DIR}/manual.pdf
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Generating documentation with pandoc"
)

