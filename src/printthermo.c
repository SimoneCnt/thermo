
/*
    Print the evaluated internal energy, entropy and free energy.

    Copyright (C) 2014, 2015 Simone Conti
    Copyright (C) 2015 Université de Strasbourg
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <thermo.h>

void 
thermo_printthermo(const Thermo *A, int onlyInt) 
{

    if (onlyInt!=1) {
    fprintf(fpout, "Extensive quantities:\n");
    fprintf(fpout, "            Elec      Trans        Rot      VibCl      VibQm      TotCl    TotQm \n");
    fprintf(fpout, "   U  %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f kcal\n",
                    A->U_elec, A->U_tr, A->U_rot, A->U_vibcl, A->U_vibqm, A->U_totcl, A->U_totqm);
    fprintf(fpout, "   S  %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f  cal\n",
                    A->S_elec, A->S_tr, A->S_rot, A->S_vibcl, A->S_vibqm, A->S_totcl, A->S_totqm);
    fprintf(fpout, " -TS  %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f kcal\n", 0.0,
                    -A->T*A->S_tr/1000.0, -A->T*A->S_rot/1000.0, -A->T*A->S_vibcl/1000.0,
                    -A->T*A->S_vibqm/1000.0, -A->T*A->S_totcl/1000.0, -A->T*A->S_totqm/1000.0);
    fprintf(fpout, "   F  %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f kcal\n",
                    A->F_elec, A->F_tr, A->F_rot, A->F_vibcl, A->F_vibqm, A->F_totcl, A->F_totqm);
    fprintf(fpout, "\n");
    }
    fprintf(fpout, "Intensive (molar) quantities:\n");
    fprintf(fpout, "            Elec      Trans        Rot      VibCl      VibQm      TotCl      TotQm \n");
    fprintf(fpout, "   Um %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f kcal/mol\n",
                    A->Um_elec, A->Um_tr, A->Um_rot, A->Um_vibcl, A->Um_vibqm, A->Um_totcl, A->Um_totqm);
    fprintf(fpout, "   Sm %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f  cal/mol\n",
                    A->Sm_elec, A->Sm_tr, A->Sm_rot, A->Sm_vibcl, A->Sm_vibqm, A->Sm_totcl, A->Sm_totqm);
    fprintf(fpout, " -TSm %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f kcal/mol\n", 0.0,
                    -A->T*A->Sm_tr/1000.0, -A->T*A->Sm_rot/1000.0, -A->T*A->Sm_vibcl/1000.0,
                    -A->T*A->Sm_vibqm/1000.0, -A->T*A->Sm_totcl/1000.0, -A->T*A->Sm_totqm/1000.0);
    fprintf(fpout, "   Fm %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f kcal/mol\n",
                    A->Fm_elec, A->Fm_tr, A->Fm_rot, A->Fm_vibcl, A->Fm_vibqm, A->Fm_totcl, A->Fm_totqm);
    fprintf(fpout, "\n");

    fprintf(fpout, "Zero point vibrational energy:  %10.3f kcal/mol\n", A->ZPE);
    
    if (A->qm_corr!=0) {
    fprintf(fpout, "Vibrational quantum correction: %10.3f kcal/mol\n", A->qm_corr);
    }

    if (!isnan(A->solvent.density)) {
        double *res = A->results;
        double dStot_omega = res[THERMO_S_EASYSOLV_TR]+res[THERMO_S_EASYSOLV_ROT]+res[THERMO_S_EASYSOLV_OMEGA];
        double dStot_epsilon = res[THERMO_S_EASYSOLV_TR]+res[THERMO_S_EASYSOLV_ROT]+res[THERMO_S_EASYSOLV_EPSILON];
        double dStot_alpha = res[THERMO_S_EASYSOLV_TR]+res[THERMO_S_EASYSOLV_ROT]+res[THERMO_S_EASYSOLV_ALPHA];
        fprintf(fpout, "\n");
        fprintf(fpout, "Solvation Entropy\n");
        fprintf(fpout, "%-20s %10s %10s %10s %10s %10s %10s             %10s %10s %10s         \n", "", "dS_trans", "dS_rot", "dS_cav", "dS_tot", "S_totcl", "S_totqm", "-TdS_tot", "F_totcl", "F_totqm");
        fprintf(fpout, "%-20s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f cal/molK    %10.3f %10.3f %10.3f kcal/mol\n", "EasySolv_Omega", res[THERMO_S_EASYSOLV_TR], res[THERMO_S_EASYSOLV_ROT],
            res[THERMO_S_EASYSOLV_OMEGA], dStot_omega, dStot_omega+A->Sm_totcl, dStot_omega+A->Sm_totqm,
            -A->T*dStot_omega/1000.0, A->Fm_totcl-A->T*dStot_omega/1000.0, A->Fm_totqm-A->T*dStot_omega/1000.0);
        fprintf(fpout, "%-20s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f cal/molK    %10.3f %10.3f %10.3f kcal/mol\n", "EasySolv_Epsilon", res[THERMO_S_EASYSOLV_TR], res[THERMO_S_EASYSOLV_ROT],
            res[THERMO_S_EASYSOLV_EPSILON], dStot_epsilon, dStot_epsilon+A->Sm_totcl, dStot_epsilon+A->Sm_totqm,
            -A->T*dStot_epsilon/1000.0, A->Fm_totcl-A->T*dStot_epsilon/1000.0, A->Fm_totqm-A->T*dStot_epsilon/1000.0);
        fprintf(fpout, "%-20s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f cal/molK    %10.3f %10.3f %10.3f kcal/mol\n", "EasySolv_Eps,Alpha", res[THERMO_S_EASYSOLV_TR], res[THERMO_S_EASYSOLV_ROT],
            res[THERMO_S_EASYSOLV_ALPHA], dStot_alpha, dStot_alpha+A->Sm_totcl, dStot_alpha+A->Sm_totqm,
            -A->T*dStot_alpha/1000.0, A->Fm_totcl-A->T*dStot_alpha/1000.0, A->Fm_totqm-A->T*dStot_alpha/1000.0);
    }

    fprintf(fpout, "\n");

    return;
}

