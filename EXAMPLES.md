
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


