
/*
    Intialize to default values a Thermo structure.

    Copyright (C) 2014-2017 Simone Conti
    Copyright (C) 2015 Université de Strasbourg
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <thermo.h>

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
    A->hessfile = NULL;
    return;
}

