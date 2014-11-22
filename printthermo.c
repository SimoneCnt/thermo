/******************************************************************************
 *
 *  thermo/printthermo.c
 *  Print to stdout the evaluated internal energy, entropy and free energy
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
 * @brief Print to stdout the evaluated internal energy, entropy and free energy
 *
 * Print to stdout the evaluated internal energy, entropy and free energy
 * 
 * @param[in] A Pointer to an initialized @c Thermo structure
 * @param[in] onlyInt If set to one, only the intensive quantities are printed
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermo.h"

void 
thermo_printthermo(const Thermo *A, int onlyInt) 
{

    if (onlyInt!=1) {
    printf ("Extensive quantities:\n");
    printf("           Elec      Trans        Rot        Vib      Total \n");
    printf("  U  %10.3f %10.3f %10.3f %10.3f %10.3f  kcal\n", A->U_elec, A->U_tr, A->U_rot, A->U_vib, A->U_tot);
    printf("  S  %10.3f %10.3f %10.3f %10.3f %10.3f   cal\n", A->S_elec, A->S_tr, A->S_rot, A->S_vib, A->S_tot);
    printf("-TS  %10.3f %10.3f %10.3f %10.3f %10.3f  kcal\n", 0.0, -A->T*A->S_tr/1000.0, -A->T*A->S_rot/1000.0, -A->T*A->S_vib/1000.0, -A->T*A->S_tot/1000.0);
    printf("  F  %10.3f %10.3f %10.3f %10.3f %10.3f  kcal\n", A->F_elec, A->F_tr, A->F_rot, A->F_vib, A->F_tot);
    printf("\n");
    }
    printf ("Intensive (molar) quantities:\n");
    printf("           Elec      Trans        Rot        Vib      Total \n");
    printf("  Um %10.3f %10.3f %10.3f %10.3f %10.3f  kcal/mol\n", A->Um_elec, A->Um_tr, A->Um_rot, A->Um_vib, A->Um_tot);
    printf("  Sm %10.3f %10.3f %10.3f %10.3f %10.3f   cal/mol\n", A->Sm_elec, A->Sm_tr, A->Sm_rot, A->Sm_vib, A->Sm_tot);
    printf("-TSm %10.3f %10.3f %10.3f %10.3f %10.3f  kcal/mol\n", 0.0, -A->T*A->Sm_tr/1000.0, -A->T*A->Sm_rot/1000.0, -A->T*A->Sm_vib/1000.0, -A->T*A->Sm_tot/1000.0);
    printf("  Fm %10.3f %10.3f %10.3f %10.3f %10.3f  kcal/mol\n", A->Fm_elec, A->Fm_tr, A->Fm_rot, A->Fm_vib, A->Fm_tot);
    printf("\n");

    return;
}

