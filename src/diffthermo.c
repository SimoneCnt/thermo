/******************************************************************************
 *
 *  thermo/diffthermo.c
 *  Given two systems and the stechiometric coefficents calculates the 
 *  thermodynamic quantities for the reaction aA<->bB
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
 * @brief Given two systems and the stechiometric coefficents calculates the 
 *        thermodynamic quantities for the reaction aA<->bB
 *
 * Given two systems and the stechiometric coefficents calculates the 
 * thermodynamic quantities for the reaction aA<->bB. Evaluates D = nB*B - nA*A
 * 
 * @param[in] A Pointer to an initialized @c Thermo structure
 * @param[in] B Pointer to an initialized @c Thermo structure
 * @param[in] nA Stechiometric coefficient for structure @c A
 * @param[in] nB Stechiometric coefficient for structure @c B
 * @param[out] D Pointer to an initialized @c Thermo structure
 *
 */

#include <stdio.h>
#include "thermo.h"

void
thermo_diffthermo(const Thermo *A, const Thermo *B, int nA, int nB, Thermo *D)
{
    D->q_elec = 0; // nB * B->q_elec  - nA * A->q_elec    ; 
    D->q_tr   = 0; // nB * B->q_tr    - nA * A->q_tr      ;
    D->q_rot  = 0; // nB * B->q_rot   - nA * A->q_rot     ;
    D->q_vibcl= 0; // nB * B->q_vib   - nA * A->q_vib     ;
    D->q_vibqm= 0; // nB * B->q_vib   - nA * A->q_vib     ;
    D->q_totcl= 0; // nB * B->q_tot   - nA * A->q_tot     ;
    D->q_totqm= 0; // nB * B->q_tot   - nA * A->q_tot     ;

    D->S_elec = 0; // nB * B->S_elec  - nA * A->S_elec    ;
    D->S_tr   = 0; // nB * B->S_tr    - nA * A->S_tr      ;
    D->S_rot  = 0; // nB * B->S_rot   - nA * A->S_rot     ;
    D->S_vibcl= 0; // nB * B->S_vib   - nA * A->S_vib     ;
    D->S_vibqm= 0; // nB * B->S_vib   - nA * A->S_vib     ;
    D->S_totcl= 0; // nB * B->S_tot   - nA * A->S_tot     ;
    D->S_totqm= 0; // nB * B->S_tot   - nA * A->S_tot     ;

    D->U_elec = 0; // nB * B->U_elec  - nA * A->U_elec    ;
    D->U_tr   = 0; // nB * B->U_tr    - nA * A->U_tr      ;
    D->U_rot  = 0; // nB * B->U_rot   - nA * A->U_rot     ;
    D->U_vibcl= 0; // nB * B->U_vib   - nA * A->U_vib     ;
    D->U_vibqm= 0; // nB * B->U_vib   - nA * A->U_vib     ;
    D->U_totcl= 0; // nB * B->U_tot   - nA * A->U_tot     ;
    D->U_totqm= 0; // nB * B->U_tot   - nA * A->U_tot     ;

    D->F_elec = 0; // nB * B->F_elec  - nA * A->F_elec    ;
    D->F_tr   = 0; // nB * B->F_tr    - nA * A->F_tr      ;
    D->F_rot  = 0; // nB * B->F_rot   - nA * A->F_rot     ;
    D->F_vibcl= 0; // nB * B->F_vib   - nA * A->F_vib     ;
    D->F_vibqm= 0; // nB * B->F_vib   - nA * A->F_vib     ;
    D->F_totcl= 0; // nB * B->F_tot   - nA * A->F_tot     ;
    D->F_totqm= 0; // nB * B->F_tot   - nA * A->F_tot     ;

    D->Sm_elec  = nB * B->Sm_elec  - nA * A->Sm_elec  ;
    D->Sm_tr    = nB * B->Sm_tr    - nA * A->Sm_tr    ;
    D->Sm_rot   = nB * B->Sm_rot   - nA * A->Sm_rot   ;
    D->Sm_vibcl = nB * B->Sm_vibcl - nA * A->Sm_vibcl ;
    D->Sm_vibqm = nB * B->Sm_vibqm - nA * A->Sm_vibqm ;
    D->Sm_totcl = nB * B->Sm_totcl - nA * A->Sm_totcl ;
    D->Sm_totqm = nB * B->Sm_totqm - nA * A->Sm_totqm ;

    D->Um_elec  = nB * B->Um_elec  - nA * A->Um_elec  ;
    D->Um_tr    = nB * B->Um_tr    - nA * A->Um_tr    ;
    D->Um_rot   = nB * B->Um_rot   - nA * A->Um_rot   ;
    D->Um_vibcl = nB * B->Um_vibcl - nA * A->Um_vibcl ;
    D->Um_vibqm = nB * B->Um_vibqm - nA * A->Um_vibqm ;
    D->Um_totcl = nB * B->Um_totcl - nA * A->Um_totcl ;
    D->Um_totqm = nB * B->Um_totqm - nA * A->Um_totqm ;

    D->Fm_elec  = nB * B->Fm_elec  - nA * A->Fm_elec  ;
    D->Fm_tr    = nB * B->Fm_tr    - nA * A->Fm_tr    ;
    D->Fm_rot   = nB * B->Fm_rot   - nA * A->Fm_rot   ;
    D->Fm_vibcl = nB * B->Fm_vibcl - nA * A->Fm_vibcl ;
    D->Fm_vibqm = nB * B->Fm_vibqm - nA * A->Fm_vibqm ;
    D->Fm_totcl = nB * B->Fm_totcl - nA * A->Fm_totcl ;
    D->Fm_totqm = nB * B->Fm_totqm - nA * A->Fm_totqm ;

    D->ZPE      = nB * B->ZPE      - nA * A->ZPE      ;

    /* Vibrational quantum correction ~ see M.Cecchini, JCTC 2015 */
    if (nA==1 && nB==1) {
        D->qm_corr = D->Fm_totqm - D->Fm_totcl;
    } else {
        D->qm_corr = 0;
    }

    return;
}

