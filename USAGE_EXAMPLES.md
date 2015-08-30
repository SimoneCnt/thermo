
Usage Examples
==============

Here a small description on how to use the code is given. Inside the `examples` 
directory it is possible to find some ready to use input files for some more 
applications. A reference output is normally present with the `.ref` suffix.


Water
-----

The easiest system we can take as example is water. Here a sample input script
to evaluate the thermodynamic properties of an hypotetical box of 22.465 liters 
which contains 1 mol of gas water at 298.15 kelvin (1 atm pressure). Each water 
molecule will have three translational degrees of freedom, three rotational 
(with the associated moments of inertia), and \f$3N-6=3\f$ vibrations (and 
associated frequencies). The simmetry number is also reported: since water 
symmetry pointgroup is C2v, the symmetry number is 2. The mass is also specified. 

    # Set temperature in Kelvin
    T 298.15

    # Set number of mols and volume (1 atm pressure)
    n 1
    V 22.465

    # Mass in g/mol
    m 18.01528

    # Translational degrees of freedom
    t 3

    # Rotational degrees of freedom and moments of inertia in g/mol Ang^2
    r 3
    1.7704
    0.6169
    1.1535

    # Symmetry number
    s 2

    # Number of vibrations and their frequencies in cm-1
    v 3
    1635.618
    3849.420
    3974.869

Save this input as e.g. `water.inp` and run it with `thermo -A water.inp` 
The output will print all parsed options and at the end all thermodynamic 
quantites, like energy, entropy, free energy, and chemical potential, are 
printed. For example, the evaluated translational molar entropy is evaluated 
to be 32.452 cal mol<sup>-1</sup> K<sup>-1</sup> and the zero point vibrational
energy to be 13.524 kcal mol<sup>-1</sup>.

Input and reference output can be found in the `water` subdirectory.


Dimerization of Insulin
-----------------------

A second example in based on the reference paper by B. Tidor and 
M. Karplus \cite tidor1994contribution, published in 1994:

> B. Tidor and M.Karplus, "The Contribution of Vibrational Entropy to 
> Molecular Association: The Dimerization of Insulin", Journal of 
> Molecular Biology, 238(3) 1994 pp. 405-414

In this work the effect of vibrations in the stabilization of the insulin 
dimer is discussed as example of how the vibrational entropy can partially 
counterbalance the net lost in translational and rotational degrees of freedom.
Translational and rotational contribution are calculated in good agreement 
with the data reported in the paper, while the vibrational ones are incorrect 
due to the limited number of frequences reported in the main text.

The input files are inside the `insulin` subdirectory. Since we want to study
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

 -  alanine dipeptide (`diala` directory) in its c7 equatorial (c7eq) and c7 
    axial (c7ax) conformations;
 -  beta hairpin from protein G (`bhp`) in the native beta-harpin (bhp1)
    and in a three-stranded beta sheet (bhp2) conformations;
 -  the converter of the biomolecular motor myosin VI (`conv`) in the 
    pre-powerstroke (pps) and rigor-like (rig) conformations.

\latexonly
\begin{figure}[ht]
\centering
\includegraphics[width=0.5\textwidth]{../mc_qmcorr.pdf}
\caption{Representation of the three examples to study the quantum correction 
to the conformational free energy difference: the alanine dipeptide, the beta 
hairpin of protein G and the converter of myosin VI. Reproduced 
from~\cite{cecchini2015quantum}.}
\label{fig_mc_qmcorr}
\end{figure}
\endlatexonly

These examples, see Fig. \latexonly\ref{fig_mc_qmcorr}\endlatexonly, are taken 
form the work of M. Cecchini, JCTC 2015 \cite cecchini2015quantum , where they 
were used to develop a quantum correction to the classical conformational free 
energy difference. 

> M. Cecchini, "Quantum Corrections to the Free Energy Difference between
> Peptides and Proteins Conformers", Journal of Chemical Theory and 
> Computation, 2015 ASAP

The quantum correction is by default calculated by Thermo; for example, for 
the myosin converter

    thermo -A rig.inp -B pps.inp --stechio 1:1

Thermo calculates a quantum correction of 1.09 kcal mol<sup>-1</sup>, in 
perfect agreement with the Cecchini paper. 


