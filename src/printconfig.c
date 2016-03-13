/******************************************************************************
 *
 *  thermo/printconfig.c
 *  Print to stdout the parsed quantities from a Thermo input file.
 *
 *  Copyright (C) 2014, 2015 Simone Conti
 *  Copyright (C) 2015 Université de Strasbourg
 *
 *  This file is part of Thermo
 * 
 *  Thermo is free software: you can redistribute it and/or modify it under the
 *  terms of the GNU General Public License as published by the Free Software 
 *  Foundation, either version 3 of the License, or (at your option) any later 
 *  version.
 *
 *  Thermo is distributed in the hope that it will be useful, but WITHOUT ANY 
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 *  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
 *  details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/**
 * @ingroup modthermo
 *
 * Print to stdout the parsed quantities from a Thermo input file.
 * 
 * @param[in] A Pointer to an initialized @c Thermo structure
 *
 */


#include <stdio.h>
#include "thermo.h"

void
thermo_printconfig(const Thermo *A) 
{
    int i;

    printf("Parsed thermodynamic quantities:\n");
    printf("   Temperature [K]:           %g\n",A->T);
    printf("   Number of moles [mol]:     %g\n",A->n);
    printf("   Volume [dm^3]:             %g\n",A->V);
    printf("   Concentration [M]:         %g\n", A->n/A->V);
    printf("   Molecular mass [g/mol]:    %g\n",A->m);
    printf("   Molar energy [kcal/mol]:   %.6f\n",A->E);
    printf("   Degree of freedom:\n");
    printf("      translational:          %d\n",A->t);
    printf("      rotational:             %d\n",A->r);
    printf("         moments of inerzia [g/mol/A^2]:\n");
    for (i=0; i<A->r; i++) { 
    printf("            %.6f\n",A->I[i]); 
    }
    printf("         symmetry number:     %d\n",A->s);
    printf("      vibrational modes:      %d\n",A->v);
    printf("         frequencies [1/cm]: \n");
    for (i=0; i<A->v; i++) {
        printf("%11.6f  ",A->nu[i]);
        if ((i+1)%6==0 && i+1!=A->v) {printf("\n");}
    }
    printf("\n");
    printf("\n");
    return;
}

