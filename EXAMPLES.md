
Examples
========

Here a small description on how to use the code is given. Inside the `examples` 
directory it is possible to find some ready to use input files for some more 
applications. A reference output is normally present with the `.ref` suffix.


Water
-----

The easiest system we can take as example is water. Here a sample input script
to evaluate the thermodynamic properties of an hypotetical box of 22.465 liters 
which contains 1 mol of gas water at 298.15 kelvin (1 atm pressure). Each water 
molecule will have three translational degrees of freedom, three rotational 
(with the associated moments of inertia), and $3N-6=3$ vibrations (and
associated frequencies). The simmetry number is also reported: since water 
symmetry pointgroup is C2v, the symmetry number is 2. The mass is also specified. 

    # Set temperature in Kelvin
    temperature = 298.15

    # Set number of mols and volume (1 atm pressure)
    nmols = 1
    volume = 22.465

    # Mass in g/mol
    mass = 18.01528

    # Translational degrees of freedom
    translational = 3

    # Rotational degrees of freedom and moments of inertia in g/mol Ang^2
    rotational = 3
    1.7704
    0.6169
    1.1535

    # Symmetry number
    sigma = 2

    # Number of vibrations and their frequencies in cm-1
    vibrations = 3
    1635.618
    3849.420
    3974.869

Save this input as e.g. `water.inp` and run it with `thermo -A water.inp` 
The output will print all parsed options and at the end all thermodynamic 
quantites, like energy, entropy, free energy, and chemical potential, are 
printed. For example, the evaluated translational molar entropy is evaluated 
to be 32.452 cal mol<sup>-1</sup> K<sup>-1</sup> and the zero point vibrational
energy to be 13.524 kcal mol<sup>-1</sup>.

Input and reference output can be found in the `examples/water` subdirectory.


Dimerization of Insulin
-----------------------

A second example in based on the reference paper by B. Tidor and 
M. Karplus [@tidor1994contribution], published in 1994:

> B. Tidor and M.Karplus, "The Contribution of Vibrational Entropy to 
> Molecular Association: The Dimerization of Insulin", Journal of 
> Molecular Biology, 238(3) 1994 pp. 405-414

In this work the effect of vibrations in the stabilization of the insulin 
dimer is discussed as example of how the vibrational entropy can partially 
counterbalance the net lost in translational and rotational degrees of freedom.
Translational and rotational contribution are calculated in good agreement 
with the data reported in the paper, while the vibrational ones are incorrect 
due to the limited number of frequences reported in the paper itself.

The input files are inside the `examples/insulin` subdirectory. Since we want to study
the dimerization reaction, it is possible to run at the same time both the 
monomer and the dimer, and, thanks to the `--stechio`` command line option,
the differences for the dimerization reaction are automatically printed. 
The Thermo command line can thus be:

    thermo -A monomer.inp -B dimer.inp --stechio 2:1 > insulin.out

The reference output is `insulin.out.ref`. At the moment only reactions between
two species can be studied (so only `-A` and `-B` command are available).


Free Energy Difference between Peptides and Proteins Conformers
---------------------------------------------------------------

Three examples are included to study the difference in free energy between 
peptides and proteins conformes: 

 -  alanine dipeptide (`examples/diala` directory) in its c7 equatorial (c7eq) and c7 
    axial (c7ax) conformations;
 -  beta hairpin from protein G (`examples/bhp`) in the native beta-harpin (bhp1)
    and in a three-stranded beta sheet (bhp2) conformations;
 -  the converter of the biomolecular motor myosin VI (`examples/conv`) in the 
    pre-powerstroke (pps) and rigor-like (rig) conformations.

![Representation of the three examples to study the quantum correction to the conformational free energy difference: the alanine dipeptide, the beta hairpin of protein G and the converter of myosin VI. Reproduced from ref. @cecchini2015quantum.\label{fig_mc_qmcorr}](doc/mc_qmcorr.pdf){width=70%}

These examples, see Fig. \ref{fig_mc_qmcorr}, are taken
form the work of M. Cecchini, JCTC 2015 [@cecchini2015quantum], where they
were used to develop a quantum correction to the classical conformational free
energy difference.

The quantum correction is by default calculated by Thermo; for example, for 
the myosin converter

    thermo -A rig.inp -B pps.inp --stechio 1:1

Thermo calculates a quantum correction of 1.09 kcal mol<sup>-1</sup>, in 
perfect agreement with the Cecchini paper. 


Solvation Entropy
-----------------

Thermo implements the "Solvation Entropy Made Simple" approach of Alejandro J.
Garza to compute solvation entropies of small molecules in different solvents
[@garza2019solvation]. This approach include three different approximations of
the solvation entropy, which are called in this code `Omega`, `Epsilon` and
`Eps,Alpha`.

To compute solvation entropies you need two input files, one for the molecule
in the gas state, one for the molecule in solution. For example, for the
solvation entropy of methanol in water the input for the methanol in gas would
look like:

```
> # Methanol gas
> temperature = 298.15
> pressure = 0.9869233
> translations = 3
> rotations = 3
> 3.98199904
> 20.50103433
> 21.27739477
> mass = 32.04
> sigma = 1
> energy = 0.0
> vibrations = 0
```

For the solution:

```
> # Methanol in water
> temperature = 298.15
> concentration = 1 unit M
> translations = 3
> rotations = 3
> 3.98199904
> 20.50103433
> 21.27739477
> mass = 32.04
> sigma = 1
> energy = 0.0
> vibrations = 0
> vvdw = 35.24161
> bbox = 119.5194
> rgyr = 1.1944
> solvent = water
```

The differences are in the concentration/pressure (1bar in gas phase, 1M in
solution), and in the addition of 4 keywords for the solution: the van der
Waals volume (vvdw), the surface of the bounding box (bbox), the radius of
gyration (rgyr), and the name of the solvent. Running the first script you
would get a total molar entropy of 53.373 cal/mol/K (Sgas). Running the second
the total molar entropy becomes 46.994 cal/mol/K (Ssolv), but three more lines
are present in the output, showing the solvation Omega, Epsilon and Eps,Alpha
solvation entropies.

The columns `dS_trans` and `dS_rot` contain the corrections to the
translational and rotational entropies, the column `dS_cav` instead contains
the cavity entropy computed in the three different methods.  Follow `dS_tot`
with the total solvation entropy correction, and `S_totcl` and `S_totqm` which
sums the solvation entropy corrections to the total entropy. Last the column
`-TdS_tot`, `F_totcl`, `F_totqm` contain the free energies. For methanol in
water, the solvation entropies corrections (dSsolv) are -21.066, -20.095 and
-19.474 cal/mol/K for the three methods.

The final solvation entropy can be obtained as (Ssolv + dSsolv) - Sgas,
obtaining -27.445, -26.474, and -25.853 cal/mol/K. The experimental value is
-27.2 cal/mol/K.

### Condensation entropy and vaporization enthalpy

In a very similar way it is possible to compute entropy of moving from the gas
phase to the liquid (condensation) and the vaporization enthalpy.  For the
condensation entropy, the final state is a molecule in solution in a solvent
composed by the same molecule, at a concentration that can be derived from the
density of the liquid. For example, for the condensation entropy of gas
methanol into liquid methanol, it is enough to change the concentration in the
previous input file from `1M` to `0.796 unit g/ml`, and change the solvent line
from water to methanol.

If the condensation -- or better the evaporation -- entropy is computed at the
boiling temperature (Tb) of the liquid, we can take advantage of the fact that
at Tb the vaporization free energy is equal to zero, which implies the
vaporization enthalpy is equal to the vaporization entropy divided by the
boiling temperature. To compute the vaporization entropy at Tb, change change
the temperature field in the input files of both the gas phase and the liquid
phase.

You can run these examples from the methanol directory:

```
# Solvation entropy of methanol from 1bar to 1M in water
thermo -A methanol-gas.thermo -B methanol-water.thermo \
    --stechio 1:1 --raw -o solution-water.out

# Vaporization entropy of methanol from the liquid at 25C
thermo -A methanol-gas.thermo -B methanol-liq.thermo \
    --stechio 1:1 --raw -o vaporization.out

# Vaporization entropy of methanol from the liquid at boiling temperature
thermo -A methanol-gas-tb.thermo -B methanol-liq-tb.thermo \
    --stechio 1:1 --raw -o vaporization-tb.out.ref
```


