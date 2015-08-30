
Examples
========

Inside the `examples` directory it is possible to find some ready to use input 
files for some published applications. A reference output is present.


Dimerization of Insulin
-----------------------

This first example in based on the reference paper by B. Tidor and 
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

The input files are inside the `insulin` subdirectory. They can be run with:

    thermo -A monomer.inp -B dimer.inp --stechio 2:1 > insulin.out

The reference output is `insulin.out.ref`


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


