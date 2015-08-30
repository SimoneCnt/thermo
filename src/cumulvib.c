/******************************************************************************
 *
 *  thermo/cumulvib.c
 *  Print the vibrational free energy in a cumulative way as a function of the 
 *  vibrational frequencies
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
 * Print the vibrational free energy in a cumulative way as a function
 * of the vibrational frequencies. This function generates two file called 
 * @c fname.k.dat and @c fname.f.dat . Both contain the cumulative vibrational 
 * free energy: the first as a function of the number of modes, the second as a 
 * function of the frequency.
 * 
 * @param[in] A     Pointer to an initialized @c Thermo structure
 * @param[in] fname Base filename to save the cumulative vibrational free energy
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "thermo.h"

void 
thermo_cumulvib(const Thermo *A, const char *fname) 
{

    int i;
    double Ftot=0, F=0;
    char fname_k[128], fname_f[128];
    strcpy(fname_k, fname);
    strcat(fname_k, ".k.dat");
    strcpy(fname_f, fname);
    strcat(fname_f, ".f.dat");

    FILE *fpk = fopen(fname_k, "w");
    FILE *fpf = fopen(fname_f, "w");
    
	for (i=0; i<A->v; i++) {
        F = A->Fm_vib_cumul_cl_k[i];
        Ftot += F;
        fprintf(fpk, "%12.4f   %12.4f   %12.4f\n", A->nu[i], F, Ftot);
	}

    Ftot = 0;
    for (i=0; i<A->nu_np; i++) {
        F = A->Fm_vib_cumul_cl[i];
        Ftot += F;
        fprintf(fpf, "%lf %lf %lf\n", A->dnu*i, F, Ftot);
    }

    fclose(fpf);
    fclose(fpk);

    return;
}

