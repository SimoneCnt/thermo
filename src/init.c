/******************************************************************************
 *
 *  thermo/init.c
 *  Intialize to default values a Thermo structure
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
 * Initialize to defualt values a Thermo structure. In particular: temperature 
 * sets to 300K, number of moles to 1, volume to 1 liter, and accuracy in 
 * vibrational to 1cm<sup>-1</sup>.
 * 
 * @param[in,out] A Pointer to an initialized @c Thermo structure
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermo.h"

/* Initialize thermo structure */
void 
thermo_init(Thermo *A) 
{

    /* Set everything to zero */
    static Thermo  tmp;
    *A=tmp;
   
    /* and set some defaults */
    A->T   = 300.0;   /* Temperature */
    A->n   = 1.0;     /* Number of moles */
    A->V   = 1.0;     /* Volume */
    A->dnu = 1.0;     /* Accuracy in vibrational spectra */
    A->nu_np = (int)lrint(ceil(4000.0/A->dnu));
    A->vdos=NULL;
    return;
}

