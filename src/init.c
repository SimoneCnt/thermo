
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
    A->s   = 1;
    A->pressure = -1;   /* Pressure (set negative as flag) */
    A->dnu = 1.0;     /* Accuracy in vibrational spectra */
    A->nu_np = (int)lrint(ceil(4000.0/A->dnu));
    A->vdos=NULL;
    A->hessfile = NULL;

    A->solute_volume = NAN;
    A->rgyr_m = NAN;
    A->asa_m = NAN;

    /* To compute solvation entropy: */
    A->solvent.name = "";
    A->solvent.vvdw = NAN;
    A->solvent.mass = NAN;
    A->solvent.density = NAN;
    A->solvent.acentricity = NAN;
    A->solvent.permittivity = NAN;
    A->solvent.expansion = NAN;
    A->solvent.rgyr = NAN;
    A->solvent.bbox = NAN;

    /* Vector containing all computed thermodynamic quantities */
    A->results = NULL;


    return;
}

