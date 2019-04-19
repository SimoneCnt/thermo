
/*
    Print the parsed quantities from a Thermo input file.

    Copyright (C) 2014-2017 Simone Conti
    Copyright (C) 2015 Université de Strasbourg
*/

#include <stdio.h>
#include <math.h>
#include <thermo.h>

void
thermo_printconfig(const Thermo *A) 
{
    int i;

    fprintf(fpout, "Parsed thermodynamic quantities:\n");
    fprintf(fpout, "   Temperature [K]:           %g\n",A->T);
    fprintf(fpout, "   Number of moles [mol]:     %g\n",A->n);
    fprintf(fpout, "   Volume [dm^3]:             %g\n",A->V);
    if (A->pressure>0.0) {
        fprintf(fpout, "   Pressure [atm]:            %g\n", A->pressure);
    } else {
        fprintf(fpout, "   Concentration [M]:         %g\n", A->n/A->V);
    }
    fprintf(fpout, "   Molecular mass [g/mol]:    %g\n",A->m);
    fprintf(fpout, "   Molar energy [kcal/mol]:   %.6f\n",A->E);
    fprintf(fpout, "   Degree of freedom:\n");
    fprintf(fpout, "      translational:          %d\n",A->t);
    fprintf(fpout, "      rotational:             %d\n",A->r);
    fprintf(fpout, "         moments of inerzia [g/mol/A^2]:\n");
    for (i=0; i<A->r; i++) { 
    fprintf(fpout, "            %.6f\n",A->I[i]);
    }
    fprintf(fpout, "         symmetry number:     %d\n",A->s);
    fprintf(fpout, "      vibrational modes:      %d\n",A->v);
    fprintf(fpout, "         frequencies [1/cm]: \n");
    for (i=0; i<A->v; i++) {
        fprintf(fpout, "%11.6f  ",A->nu[i]);
        if ((i+1)%6==0 && i+1!=A->v) {fprintf(fpout, "\n");}
    }
    fprintf(fpout, "\n");

    if (!isnan(A->solute_volume))     fprintf(fpout, "   Solute vdw volume:          %g\n", A->solute_volume);
    if (!isnan(A->rgyr_m))            fprintf(fpout, "   Solute gyration radius:     %g\n", A->rgyr_m);
    if (!isnan(A->asa_m))             fprintf(fpout, "   Solute accessible surface:  %g\n", A->asa_m);

    if (!isnan(A->solvent.density)) {
        fprintf(fpout, "\n");
        thermo_print_solvent(fpout, A->solvent);
    }

    fprintf(fpout, "\n");
    return;
}

