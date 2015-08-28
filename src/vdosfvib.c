/******************************************************************************
 *
 *  thermo/delete.c
 *  Print the cumulative vibrational free energy from the vibrational density 
 *  of states
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
 * @brief Print the cumulative vibrational free energy from the vibrational 
 *        density of states
 *
 * Print the cumulative vibrational free energy from the vibrational density 
 * of states
 * 
 * @param[in] A Pointer to an initialized @c Thermo structure
 * @param[in] fname Filename where to save the cumulative vibrational free energy
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermo.h"

void 
thermo_vdosfvib(const Thermo *A, const char *fname) 
{
    int i;                          /* Counter */
    double kBT = CNS_kB * A->T;     /* kB T */
    double F, Ftot=0;
    FILE *fp = fopen(fname, "w");

	for (i=0; i<A->nu_np; i++) {
        F = -CNS_j2kcal * CNS_NA * kBT * log( kBT / ( CNS_h * (i+1) * A->dnu * CNS_C * 100.0)) * A->vdos[i];
        Ftot += F;
        fprintf(fp, "%f %f %f\n", A->dnu*i, F, Ftot);
	}

    fclose(fp);
    return;
}

