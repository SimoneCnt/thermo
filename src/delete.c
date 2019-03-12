
/*
    Delete all allocated memory inside a Thermo structure

    Copyright (C) 2014-2017 Simone Conti
    Copyright (C) 2015 Université de Strasbourg
*/

#include <stdio.h>
#include <stdlib.h>
#include <thermo.h>

void
thermo_delete(Thermo *A)
{
    if (A->nu!=NULL) {free(A->nu); A->nu=NULL;}
    if (A->I!=NULL) {free(A->I); A->I=NULL;}
    if (A->Fm_vib_cumul_cl!=NULL) {free(A->Fm_vib_cumul_cl); A->Fm_vib_cumul_cl=NULL;}
    if (A->Fm_vib_cumul_cl_k!=NULL) {free(A->Fm_vib_cumul_cl_k); A->Fm_vib_cumul_cl_k=NULL;}
    if (A->Fm_vib_cumul_qm!=NULL) {free(A->Fm_vib_cumul_qm); A->Fm_vib_cumul_qm=NULL;}
    if (A->Fm_vib_cumul_qm_k!=NULL) {free(A->Fm_vib_cumul_qm_k); A->Fm_vib_cumul_qm_k=NULL;}
    if (A->vdos!=NULL) {free(A->vdos); A->vdos=NULL;}
    if (A->hessfile) {free(A->hessfile); A->hessfile=NULL;}
    if (A->hessian) {free(A->hessian); A->hessian=NULL;}
    return;
}

