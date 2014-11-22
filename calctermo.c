/******************************************************************************
 *
 *  thermo/calcthermo.c
 *  Calculate all thermodynamic quantities based on partition function
 *
 *  Copyright (C) 2014 Simone Conti
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
 * @brief Calculate all thermodynamic quantities based on partition function
 *
 * Taken an initialized @c Thermo structure, evaluates the translational, 
 * rotational, vibrational and electronic contributions to the partition function.
 * From that, the internal energy, entropy and free energy are evaluated. 
 * 
 * @param[in,out]  A  Pointer to an initialized @c Thermo structure
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermo.h"

void 
thermo_calcthermo(Thermo *A) 
{
    int i;                          /* Counter */
    long int ndx;
    double kBT = CNS_kB * A->T;     /* kB T */
    double tmp;

    /* Logarithm of the molecular partition function */

    /* Translational */
    A->q_tr  = (A->t / 2.0 ) * log( (CNS_2PI * (A->m/(CNS_NA*1000.0)) * kBT) / CNS_h2 ) + log(A->V/1000.0);

    /* Rotational */
    A->q_rot = log(CNS_PI/(A->s*A->s)) + A->r * log( (8.0 * CNS_PI2 * kBT) / CNS_h2);
    for (i=0; i<A->r ; i++) {
        A->q_rot += log(A->I[i] * 1E-23 / CNS_NA);
    }
    A->q_rot /= 2.0;
    
    /* Vibrational */
    A->q_vib = 0;
    A->Fm_vib_cumul_cl   = malloc(A->nu_np*sizeof(double));
    A->Fm_vib_cumul_cl_k = malloc(A->v*sizeof(double));
    for (i=0; i<A->nu_np; i++) {
        A->Fm_vib_cumul_cl[i] = 0;
    }
	for (i=0; i<A->v; i++) {

        tmp = log( kBT / ( CNS_h * A->nu[i] * CNS_C * 100.0));
		A->q_vib += tmp;

        A->Fm_vib_cumul_cl_k[i] = -CNS_j2kcal * CNS_NA * kBT * tmp;

        ndx = lrint(A->nu[i]/A->dnu);
        if (ndx>=A->nu_np) {
            printf("Warning! ndx=%ld greather than nu_np=%ld\n", ndx, A->nu_np);
            ndx = A->nu_np-1;
        }
        A->Fm_vib_cumul_cl[ndx] -= CNS_j2kcal * CNS_NA * kBT * tmp;
	}

    /* Electronic */
    A->q_elec = -A->E / (CNS_NA * kBT);

    /* Total */
    A->q_tot = A->q_tr * A->q_rot * A->q_vib * A->q_elec;


    /* Free energy */
    A->F_tr    = - CNS_j2kcal * A->n * CNS_NA * kBT * (A->q_tr + 1.0 - log(A->n*CNS_NA));   /* Translational */
    A->F_rot   = - CNS_j2kcal * A->n * CNS_NA * kBT * A->q_rot;                             /* Rotational    */
    A->F_vib   = - CNS_j2kcal * A->n * CNS_NA * kBT * A->q_vib;                             /* Vibrational   */
    A->F_elec  = - A->n * CNS_NA * kBT * A->q_elec;                                         /* Electronic    */
    A->F_tot   = A->F_tr + A->F_rot + A->F_vib + A->F_elec;                                 /* Total         */

    /* Chemical potential */
    A->Fm_tr   = A->F_tr   / (A->n) + CNS_NA*CNS_j2kcal*kBT;                                /* Translational */
    A->Fm_rot  = A->F_rot  / (A->n);                                                        /* Rotational    */
    A->Fm_vib  = A->F_vib  / (A->n);                                                        /* Vibrational   */
    A->Fm_elec = A->F_elec / (A->n);                                                        /* Electronic    */
    A->Fm_tot  = A->Fm_tr + A->Fm_rot + A->Fm_vib + A->Fm_elec;                             /* Total         */

    /* Internal Energy */
    A->U_tr    = CNS_j2kcal * A->n * CNS_NA * kBT * A->t / 2.0;                             /* Translational */
    A->U_rot   = CNS_j2kcal * A->n * CNS_NA * kBT * A->r / 2.0;                             /* Rotational    */
    A->U_vib   = CNS_j2kcal * A->n * CNS_NA * kBT * A->v;                                   /* Vibrational   */
    A->U_elec  = A->n * A->E;                                                               /* Electronic    */
    A->U_tot   = A->U_tr + A->U_rot + A->U_vib + A->U_elec;                                 /* Total         */

    /* Molar internal energy */
    A->Um_tr   = A->U_tr   / (A->n);                                                        /* Translational */
    A->Um_rot  = A->U_rot  / (A->n);                                                        /* Rotational    */
    A->Um_vib  = A->U_vib  / (A->n);                                                        /* Vibrational   */
    A->Um_elec = A->U_elec / (A->n);                                                        /* Electronic    */
    A->Um_tot  = A->U_tot  / (A->n);                                                        /* Total         */

    /* Entropy */
    A->S_tr    = - 1000.0 * (A->F_tr  - A->U_tr ) / A->T;                                   /* Translational */
    A->S_rot   = - 1000.0 * (A->F_rot - A->U_rot) / A->T;                                   /* Rotational    */
    A->S_vib   = - 1000.0 * (A->F_vib - A->U_vib) / A->T;                                   /* Vibrational   */
    A->S_elec  = 0.0;                                                                       /* Electronic    */
    A->S_tot   = A->S_tr + A->S_rot + A->S_vib;                                             /* Total         */

    /* Molar entropy */
    A->Sm_tr   = A->S_tr  / A->n - 1000.0*CNS_NA*CNS_j2kcal*CNS_kB;                         /* Translational */
    A->Sm_rot  = A->S_rot / A->n;                                                           /* Rotational    */
    A->Sm_vib  = A->S_vib / A->n;                                                           /* Vibrational   */
    A->Sm_elec = 0.0;                                                                       /* Electronic    */
    A->Sm_tot  = A->Sm_tr + A->Sm_rot + A->Sm_vib;                                          /* Total         */

    return;
}

