
/*
    Print the parsed quantities from a Thermo input file.

    Copyright (C) 2014-2017 Simone Conti
    Copyright (C) 2015 Université de Strasbourg
*/

#include <stdio.h>
#include <thermo.h>

void
thermo_printconfig(const Thermo *A) 
{
    int i;

    fprintf(fpout, "Parsed thermodynamic quantities:\n");
    fprintf(fpout, "   Temperature [K]:           %g\n",A->T);
    fprintf(fpout, "   Number of moles [mol]:     %g\n",A->n);
    fprintf(fpout, "   Volume [dm^3]:             %g\n",A->V);
    fprintf(fpout, "   Concentration [M]:         %g\n", A->n/A->V);
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
    fprintf(fpout, "\n");
    return;
}

