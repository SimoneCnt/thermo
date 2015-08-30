
Install instructions
====================

The last version of the code can be downloaded from git via:

    git clone https://github.com/SimoneCnt/thermo.git thermo

which will copy the latest available version inside the `thermo` directory.
If you do not have git installed and 
 - you want to install it, take a look at <http://git-scm.com>
 - you do not want to install it, you can get the code as a zip archive directly 
    from the git project page at <https://github.com/SimoneCnt/thermo>

Once downloaded the installation is straighforward. The code is written in 
standard C and the compilation is managed by a CMake script. So to compile 
the code the following list of commands will easily do the work:

    mkdir build
    cd build
    cmake ..
    make

If you are missing CMake, you can get if from <http://www.cmake.org>. 

Once the compilation finishes, you can find the `thermo` binary inside the 
`build` directory.

This procedure has been tested on Linux (Ubuntu) and MacOsX machines. On Ubuntu, 
and probably also in other Linux distributions, you can find both CMake and 
git on the standard packaging tools. So in Ubuntu this line should smoothly 
install both applications:

    sudo apt-get install git cmake

No test has been performed on Windows machines.

