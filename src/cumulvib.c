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
#include <thermo.h>

void 
thermo_cumulvib(const Thermo *A, const char *fname) 
{

    int i;
    double Ftot_CL, F_CL;
    double Ftot_QM, F_QM;
    char *fname_k, *fname_f;

    /* Set names for the output files */
    if (fname==NULL) {
        fname = "cumul";
    }
    size_t l = strlen(fname)+7;
    fname_k = malloc(l*sizeof(char));
    fname_f = malloc(l*sizeof(char));
    if (fname_k==NULL || fname_f==NULL) {
        fprintf(stderr, "ERROR! Memory allocation failed in cumulative free energy!\n\n");
        return;
    }
    strcpy(fname_k, fname);
    strcat(fname_k, ".k.dat");
    strcpy(fname_f, fname);
    strcat(fname_f, ".f.dat");

    /* Open output files */
    FILE *fpk = fopen(fname_k, "w");
    FILE *fpf = fopen(fname_f, "w");
    if (fpk==NULL || fpf==NULL) {
        fprintf(stderr, "ERROR opening %s or %s for comulative vibrational free energy!\n\n", fname_k, fname_f);
        free(fname_k);
        free(fname_f);
        return;
    }

    /* Print the cumulative classical and quantum vibrational free energy per mode number */
    fprintf(fpk, "#k                   FvibCL         FvibQM         Delta \n");
    Ftot_CL = 0;
    Ftot_QM = 0;
	for (i=0; i<A->v; i++) {
        F_CL = A->Fm_vib_cumul_cl_k[i];
        F_QM = A->Fm_vib_cumul_qm_k[i];
        Ftot_CL += F_CL;
        Ftot_QM += F_QM;
        fprintf(fpk, "%12.4f   %12.4f   %12.4f   %12.4f \n", A->nu[i], Ftot_CL, Ftot_QM, Ftot_QM-Ftot_CL);
	}

    /* Print the cumulative classical and quantum vibrational free energy per mode frequency */
    fprintf(fpf, "#freq                FvibCL         FvibQM         Delta \n");
    Ftot_CL = 0;
    Ftot_QM = 0;
    for (i=0; i<A->nu_np; i++) {
        F_CL = A->Fm_vib_cumul_cl[i];
        F_QM = A->Fm_vib_cumul_qm[i];
        Ftot_CL += F_CL;
        Ftot_QM += F_QM;
        fprintf(fpf, "%12.4f   %12.4f   %12.4f   %12.4f \n", A->dnu*i, Ftot_CL, Ftot_QM, Ftot_QM-Ftot_CL);
    }

    /* Clean memory and return */
    fclose(fpf);
    fclose(fpk);
    free(fname_k);
    free(fname_f);
    return;
}

