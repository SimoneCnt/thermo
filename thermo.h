/******************************************************************************
 *
 *  thermo/thermo.h
 *  Intialize to default values a Thermo structure
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

#ifndef _THERMO_
#define _THERMO_

/* Constants */
#define CNS_PI     3.14159265358979323846  /* pi []                              */
#define CNS_2PI    (2.0*CNS_PI)            /* 2pi []                             */
#define CNS_PI2    (CNS_PI*CNS_PI)         /* pi^2 []                            */
#define CNS_kB     1.3806488E-23           /* Boltzmann constant [ J/K ]         */
#define CNS_h      6.62606957E-34          /* Planck constant [ Js ]             */
#define CNS_h2     (CNS_h*CNS_h)           /* Planck constant squared [ (Js)^2 ] */
#define CNS_NA     6.02214129E+23          /* Avogadro number []                 */
#define CNS_C      299792458.0             /* Speed of light [ m/s ]             */
#define CNS_cal    4.1867999409            /* One calorie in joule [ J/cal ]     */
#define CNS_j2kcal (1.0/(CNS_cal*1000.0))  /* Convert joule to kcal [ kcal/J ]   */


/* Structure which contain all input informations about a system */
typedef struct {
    double T;   /* Temperature in kelvin */
    double V;   /* Volume */
    double n;   /* Number of mols  */
    double t;	/* Number of translation degree of freedom */
    double r;	/* Number of rotational degree of freedom */
    double v;	/* Number of vibrational degree of freedom */
    double s;	/* Symmetry number */
    double m;	/* Mass of the system in g/mol */
    double E;	/* Energy of the system in kcal/mol */
    double *I;	/* Moments of inertia in g/mol/A^2 */
    double *nu;	/* Vibrational normal modes in cm-1 */
    double dnu; /* Accuracy in vibrational spectra for cumulative and vdos calculations */
    long int nu_np;  /* TODO */
    double  q_elec,  q_tr,  q_rot,  q_vib,  q_tot; /* Natural logarithm of the molecular partition function: ln(q) */
    double  S_elec,  S_tr,  S_rot,  S_vib,  S_tot; /* Entropy */
    double  U_elec,  U_tr,  U_rot,  U_vib,  U_tot; /* Internal Energy */
    double  F_elec,  F_tr,  F_rot,  F_vib,  F_tot; /* Free energy */
    double Sm_elec, Sm_tr, Sm_rot, Sm_vib, Sm_tot; /* Molar Entropy */
    double Um_elec, Um_tr, Um_rot, Um_vib, Um_tot; /* Molar Internal Energy */
    double Fm_elec, Fm_tr, Fm_rot, Fm_vib, Fm_tot; /* Molar Free energy (chemical potential) */
    double *Fm_vib_cumul_cl, *Fm_vib_cumul_qm, *Fm_vib_cumul_delta; /* Cumulative vibrational free energy per frequency */
    double *Fm_vib_cumul_cl_k, *Fm_vib_cumul_qm_k, *Fm_vib_cumul_delta_k; /* Cumulative vibrational free energy per mode */
    double *vdos; /* Vibrational density of states */
} Thermo;


void thermo_calcthermo(Thermo *A);
void thermo_cumulvib(const Thermo *A, const char *filename);
void thermo_delete(Thermo *A);
void thermo_diffthermo(const Thermo *A, const Thermo *B, int nA, int nB, Thermo *D);
void thermo_init(Thermo *A);
void thermo_printconfig(const Thermo *A);
void thermo_printthermo(const Thermo *A, int onlyInt); 
int  thermo_readthermo(Thermo *A, const char *fname);
void thermo_vdos(Thermo *A, const char *fname);
void thermo_vdosfvib(const Thermo *A, const char *fname);

#endif
