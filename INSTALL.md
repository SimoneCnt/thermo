
Install instructions
====================

The last version of the code can be downloaded from git via:

    git clone https://github.com/SimoneCnt/thermo.git thermo

which will copy the latest available version inside the `thermo` directory.  If
you do not have git installed and you want to install it, take a look at
<http://git-scm.com> Otherwise you can get the full code as a zip archive
directly from the git project page at <https://github.com/SimoneCnt/thermo>

Once downloaded the installation is straighforward. The code is written in 
standard C and the compilation is managed by a CMake script. So to compile 
the code the following list of commands will easily do the work:

    cd thermo
    mkdir build
    cd build
    cmake ..
    make release

After these commands, you should find the `thermo` binary inside the `build`
directory. To check if the compilation went file, you can issue the `make
check` command (always from inside the build directory). 

One option of Thermo is to give as input a hessian matrix to evaluate the normal
mode frequences. To activate this option you need to link Thermo to a linear
algebra (lapack) library. On MacOS the linking should be automatic (the
developer tools in MacOS already contains an efficient lapack library).  On
Linux you need first to install an efficient lapack library (e.g. MKL from
Intel, freeware), and second to specify the include directory and the link
options to the Thermo build script. You have two options: hope it will work
automatically (unlikely), or explicitly specifify the include and link
directive through the `LAPACK_INC` and `LAPACK_LINK` environmental variables.
For example, on my linux machine to link to MKL I have:

    LAPACK_LINK=-L/opt/compilers_and_libraries_2017.0.098/linux/mkl/lib/intel64 
        -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp 
        -lpthread -lm -ldl
    LAPACK_INCL=-m64 -I/opt/compilers_and_libraries_2017.0.098/linux/mkl/include

This is higly depenent on your system and which lapack library you use.

For any problem, feel free to refer to <https://github.com/SimoneCnt/thermo/issues> 
or directly to <simonecnt@gmail.com>.

