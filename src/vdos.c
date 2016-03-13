/******************************************************************************
 *
 *  thermo/vdos.c
 *  Evaluate the vibrational density of states VDOS starting from the 
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
 * Evaluate the vibrational density of states (VDOS) starting from the 
 * vibrational frequencies. The free energy is vibrational free energy is also
 * evaluated as integral over the VDOS (see Theory). In the output file, the VDOS, 
 * the free energy at each point, and the cumulative free energy is printes as a 
 * function of the frequency. To change the resolutin of the calculated VDOS, 
 * use the -n, --dnu command line option.
 *
 * @param[in] A Pointer to an initialized @c Thermo structure
 * @param[in] fname Filename where to print the VDOS and the free energy.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermo.h"

void 
thermo_vdos(Thermo *A, const char *fname)
{

    int i, ii;
    int period = 5;                 /* Period for moving average of VDOS */
    int weig, totweig;
    double accu, F, Ftot;
    double kBT = CNS_kB * A->T;
    double *vdos_tmp;
    FILE *out = fopen(fname, "w");

    /* Allocate the vdos and initialize to zero */
    A->vdos  = malloc((size_t)(A->nu_np)*sizeof(double));
    vdos_tmp = malloc((size_t)(A->nu_np)*sizeof(double));
    for (i=0; i<A->nu_np; i++) {
        A->vdos[i]  = 0;
        vdos_tmp[i] = 0;
    }

    /* Populate the vdos */
    for (i=0; i<A->v; i++) {
        ii = (int)lrint(A->nu[i]/A->dnu);
        if (ii>=A->nu_np) ii=A->nu_np-1;
        A->vdos[ii]++;
    }

    /* Do moving average */
    for (i=0; i<A->nu_np; i++) {
        totweig = 0;
        for (ii=-period; ii<=period; ii++) {
            if (i+ii>0 && i+ii<A->nu_np) {
                weig = period+1-abs(ii);
                vdos_tmp[i] += A->vdos[i+ii] * weig;
                totweig  += weig;
            }
        }
        vdos_tmp[i] /= totweig;
    }

    /* Normalize it */
    accu=0;
    for (i=0; i<A->nu_np; i++) {
        accu += vdos_tmp[i];
    }
    for (i=0; i<A->nu_np; i++) {
        A->vdos[i] = vdos_tmp[i] * A->v / accu;
    }

    /* Print the vdos, the free energy gain at each point, and the cumulative free energy */
    fprintf(out, "#freq      vdos           F        Ftot\n");
    Ftot=0;
    for (i=0; i<A->nu_np; i++) {
        F = -CNS_j2kcal * CNS_NA * kBT * log( kBT / ( CNS_h * (i+1) * A->dnu * CNS_C * 100.0)) * A->vdos[i];
        Ftot += F;
        fprintf(out, "%7.2f %+11.4e %+10.4f %+10.4f \n", A->dnu*i, A->vdos[i], F, Ftot);
    }

    free(vdos_tmp);
    fclose(out);

    return;
}
