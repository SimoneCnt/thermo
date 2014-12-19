/******************************************************************************
 *
 *  thermo/delete.c
 *  Delete all allocated memory inside a Thermo structure
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
 * @brief Delete all allocated memory inside a Thermo structure
 *
 * Delete all allocated memory inside a Thermo structure
 * 
 * @param[in,out] A Pointer to an initialized @c Thermo structure
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include "thermo.h"

void
thermo_delete(Thermo *A)
{
    if (A->nu!=NULL) {free(A->nu); A->nu=NULL;}
    if (A->I!=NULL) {free(A->I); A->I=NULL;}
    if (A->Fm_vib_cumul_cl!=NULL) {free(A->Fm_vib_cumul_cl); A->Fm_vib_cumul_cl=NULL;}
    if (A->Fm_vib_cumul_cl_k!=NULL) {free(A->Fm_vib_cumul_cl_k); A->Fm_vib_cumul_cl_k=NULL;}
    if (A->vdos!=NULL) {free(A->vdos); A->vdos=NULL;}
    return;
}

