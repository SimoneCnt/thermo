
Usage
=====

To use Thermo you need to set some command line options and to write a thermo
input file. For each molecule you wanto to study, you need to set up one thermo
input file, which contains the molecular properties used to evaluate the
entropy, the free energy and the internal energy.  The thermo input is
described first, the command line options follow.


Thermo input
------------

The thermo input file is a simple text file which constains `key = value`
pairs. Based on the "key" set, different contributions to the total free energy
are evaluated. Generally speaking (see the Theory section for more details),
the free energy (as the entropy and the internal energy) is split in an
energetic, a translational, rotational and vibrational contributions. Usually,
a molecule in solution has three rotational and three translational degrees of
freedom and 3N-6 vibrations (N is the number of atoms). This is not valid for
linear molecules or for systems with constrained degrees of freedom, like
molecules adsorbed on a surface. The keys described below allow you to define
the repartition of the 3N degrees of freedom and the molecular properties, like
mass and moment of inertia, necessary to evaluate the free energy.

The keys you can set in the thermo input file are:

* **temperature**: The temperature in kelvin. Default 300 K.
* **nmoles**: The number of moles in mol. Default 1 mol.
* **volume**: Volume in liters. Together with the number of moles defines the 
    concentration. Default 1 l.
* **energy**: The ground state energy in kcal/mol. Default 0 kcal/mol.
* **mass**: Molecular mass in g/mol. Default zero.
* **translations**: Number of translational degrees of freedom. Default 3.
* **rotations**: Number of rotational degrees of freedom. Default 3. 
    This variable has to be followed by the moments of inertia one per line, 
    in equal number of the specified number of rotational degrees of freedom, 
    expressed in g/mol/A^2^. 
* **sigma**: Symmetry number. Default 1.
* **vibrations**: Number of vibrational degrees of freedom. Default zero. 
    This variable has to be followed by the vibrational frequencies, one per line, 
    in equal number of the specified umber of vibrational degrees of freedom, 
    expressed in cm^-1^.
* **hessian**: Name of the file where a hessian matrix can be read. The matrix 
    will be internally diagonalized (if Thermo is compiled with lapack support), 
    and the normal mode frequences calculated. The hessian and vibrations 
    keywords are mutually exclusive. Default none.

For example, if you want to specify the temperature you can write in the input file:

    temperature = 310

or for the moments of inertia:

    rotations = 3
    moment1
    moment2
    moment3

The keywords are identified based only on the first four charachters, so these
three statements are actually identical:

    temperature = 300 
    temp = 300
    temperatura = 300


Command line options
--------------------

The mode of working of Thermo id defined by few command line options. The basic
is via the `--A` option followed by the name of the thermo input file. Thus
running

    thermo --A myinput.thermo 

will calculate all internal energy, entropy and free energy based on the
options set in the `myinput.thermo` imput file. The `--B` option is essetianlly
the same as `--A`. You can use both `--A` and `--B` in the same run. This is
useful in particular with the `--stechio a:b` option. The parameters `a` and
`b` are the stechiometric coefficient for the reaction $aA \rightleftharpoons
bB$. Thus via this command line

    thermo --A monomer.thermo --B dimer.thermo --stechio 2:1

thermo will evaluate the free energy of the monomer, the dimer, and the free
energy difference for the dimerization reaction.

More classical command line options, `--out outfile.out` redirect the thermo output to the `outfile.out` file, `--help` print an hopefully useful help, and `--version` print the current version of the thermo code.

Still to document: `--cumul`, `--vdos`, `--dnu`. These essentially create and write to file the vibrational density of states (VDOS) and the cumulative vibrational free energy.

