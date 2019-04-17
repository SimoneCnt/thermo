
/*
    Calculate all thermodynamic quantities based on partition function.
    Taken an initialized Thermo structure, evaluates the translational, 
    rotational, vibrational and electronic contributions to the partition function.
    From that, the internal energy, entropy and free energy are evaluated. 

    Copyright (C) 2014-2019 Simone Conti
    Copyright (C) 2015 Université de Strasbourg
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <thermo.h>

void 
thermo_calcthermo(Thermo *A) 
{
    int i;                          /* Counter */
    int ndx;
    double kBT = CNS_kB * A->T;     /* kB T */
    double tmp;

    A->results = thermo_compute(A->T, A->E, A->t, A->m, A->V, A->n, A->r, A->I, A->s, A->v, A->nu,
        A->solute_volume, A->solvent_volume, A->solvent_mass, A->density,
        A->acentricity, A->permittivity, A->thermal_expansion,
        A->rgyr_m, A->rgyr_s, A->asa_m, A->asa_s);
    double *res = A->results;

    if (!res) {
        fprintf(stderr, "thermo_calcthermo: thermo computation failed!\n");
        return;
    }

    A->q_elec  = res[THERMO_LNQ_ELEC];
    A->Fm_elec = res[THERMO_F_ELEC];
    A->Um_elec = res[THERMO_U_ELEC];
    A->Sm_elec = res[THERMO_S_ELEC];

    A->q_tr  = res[THERMO_LNQ_TR];
    A->Fm_tr = res[THERMO_F_TR];
    A->Um_tr = res[THERMO_U_TR];
    A->Sm_tr = res[THERMO_S_TR];
    
    A->q_rot  = res[THERMO_LNQ_ROT];
    A->Fm_rot = res[THERMO_F_ROT];
    A->Um_rot = res[THERMO_U_ROT];
    A->Sm_rot = res[THERMO_S_ROT];

    A->q_vibcl  = res[THERMO_LNQ_VIBCL];
    A->Fm_vibcl = res[THERMO_F_VIBCL];
    A->Um_vibcl = res[THERMO_U_VIBCL];
    A->Sm_vibcl = res[THERMO_S_VIBCL];

    A->q_vibqm  = res[THERMO_LNQ_VIBQM];
    A->Fm_vibqm = res[THERMO_F_VIBQM];
    A->Um_vibqm = res[THERMO_U_VIBQM];
    A->Sm_vibqm = res[THERMO_S_VIBQM];
    A->ZPE = res[THERMO_ZPE];

    /* Cumulative Vibrational free energy */
    A->Fm_vib_cumul_cl   = malloc((size_t)(A->nu_np)*sizeof(double));
    A->Fm_vib_cumul_cl_k = malloc((size_t)(A->v)*sizeof(double));
    A->Fm_vib_cumul_qm   = malloc((size_t)(A->nu_np)*sizeof(double));
    A->Fm_vib_cumul_qm_k = malloc((size_t)(A->v)*sizeof(double));
    for (i=0; i<A->nu_np; i++) {
        A->Fm_vib_cumul_cl[i] = 0;
        A->Fm_vib_cumul_qm[i] = 0;
    }
	for (i=0; i<A->v; i++) {
        ndx = (int)lrint(A->nu[i]/A->dnu);
        if (ndx>=A->nu_np) {
            fprintf(fpout, "Warning! ndx=%d greather than nu_np=%d\n", ndx, A->nu_np);
            ndx = A->nu_np-1;
        }
        /* Classic */
        tmp = log( kBT / ( CNS_h * A->nu[i] * CNS_C * 100.0));
        A->Fm_vib_cumul_cl_k[i] = -CNS_j2kcal * CNS_NA * kBT * tmp;
        A->Fm_vib_cumul_cl[ndx] -= CNS_j2kcal * CNS_NA * kBT * tmp;
        /* Quantic */
        tmp = ( CNS_h * A->nu[i] * CNS_C * 100.0) / (2.0*kBT);
        A->Fm_vib_cumul_qm_k[i] = -CNS_j2kcal * CNS_NA * kBT * (-log(2.0*sinh(tmp)));
        A->Fm_vib_cumul_qm[ndx] -= CNS_j2kcal * CNS_NA * kBT * (-log(2.0*sinh(tmp)));
	}

    /* Vibrational quantum correction ~ see M. Cecchini, JCTC 2015 */
    A->qm_corr = 0;

    /* Total */
    A->q_totcl = A->q_tr + A->q_rot + A->q_vibcl + A->q_elec;
    A->q_totqm = A->q_tr + A->q_rot + A->q_vibqm + A->q_elec;

    A->Fm_totcl = A->Fm_tr + A->Fm_rot + A->Fm_vibcl + A->Fm_elec;
    A->Fm_totqm = A->Fm_tr + A->Fm_rot + A->Fm_vibqm + A->Fm_elec;
    A->Um_totcl = A->Um_tr + A->Um_rot + A->Um_vibcl + A->Um_elec;
    A->Um_totqm = A->Um_tr + A->Um_rot + A->Um_vibqm + A->Um_elec;
    A->Sm_totcl = A->Sm_tr + A->Sm_rot + A->Sm_vibcl + A->Sm_elec;
    A->Sm_totqm = A->Sm_tr + A->Sm_rot + A->Sm_vibqm + A->Sm_elec;

    /* Molar quantities computed above, now doing extensive ones */

    A->F_tr    = A->n * (A->Fm_tr - CNS_NA*CNS_j2kcal*CNS_kB*A->T);
    A->F_rot   = A->n * A->Fm_rot;
    A->F_vibcl = A->n * A->Fm_vibcl;
    A->F_vibqm = A->n * A->Fm_vibqm;
    A->F_elec  = A->n * A->Fm_elec;
    A->F_totcl = A->F_tr + A->F_rot + A->F_vibcl + A->F_elec;
    A->F_totqm = A->F_tr + A->F_rot + A->F_vibqm + A->F_elec;

    A->U_tr    = A->n * A->Um_tr;
    A->U_rot   = A->n * A->Um_rot;
    A->U_vibcl = A->n * A->Um_vibcl;
    A->U_vibqm = A->n * A->Um_vibqm;
    A->U_elec  = A->n * A->Um_elec;
    A->U_totcl = A->U_tr + A->U_rot + A->U_vibcl + A->U_elec;
    A->U_totqm = A->U_tr + A->U_rot + A->U_vibqm + A->U_elec;

    A->S_tr    = A->n * (A->Sm_tr + 1000.0*CNS_NA*CNS_j2kcal*CNS_kB);
    A->S_rot   = A->n * A->Sm_rot;
    A->S_vibcl = A->n * A->Sm_vibcl;
    A->S_vibqm = A->n * A->Sm_vibqm;
    A->S_elec  = A->n * A->Sm_elec;
    A->S_totcl = A->S_tr + A->S_rot + A->S_vibcl + A->S_elec;
    A->S_totqm = A->S_tr + A->S_rot + A->S_vibqm + A->S_elec;

    return;
}

