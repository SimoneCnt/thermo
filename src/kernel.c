
/*
    Calculate thermodynamic quantities based on partition function.

    Copyright (C) 2019 Simone Conti
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <thermo.h>


#define PLANCK      6.62606957E-34              /* Planck constant [ Js ] */
#define AVOGADRO    6.02214129E+23              /* Avogadro number [mol-1] */
#define BOLTZMANN   1.3806488E-23               /* Boltzmann constant [ J/K ] */
#define LIGHTSPEED  299792458.0                 /* Speed of light [ m/s ]             */
#define J2KCALMOL   (AVOGADRO/(4.184*1000.0))   /* Convert joule to kcal/mol [ kcal/J ] */
#define PI          M_PI                        /* 3.14... */


/*
    Convert frequency in cm-1 to kelvin.
*/
double thermo_cm2kelvin(double cm) {
    return cm/(BOLTZMANN/(PLANCK*LIGHTSPEED*100.0));

}

/*
    Convert vibrational frequencies from kelvin to cm-1.
*/
double thermo_kelvin2cm(double kelvin) {
    return kelvin*(BOLTZMANN/(PLANCK*LIGHTSPEED*100.0));
}

/*
    Convert a rotational temperature to a moment of inertia in (g/mol)*A^2.
*/
double thermo_kelvin2inertia(double kelvin) {
    return ((PLANCK*PLANCK*AVOGADRO)/(8.0*PI*PI*BOLTZMANN*1E-23))/kelvin;
}

/*
    Convert rotational frequencies [cm-1] to moment of inertia [(g/mol)*A^2]
*/
double thermo_freq2inertia(double freq) {
    return thermo_kelvin2inertia(thermo_cm2kelvin(freq));
}



/*
    Ideal gas translations
    ----------------------

    Compute thermodynamic quantities for translational degrees of freedom, for an ideal gas.
    The logarithm of the partition function is:

    log(q_tr) = log [ (2 pi m kB T / h2 )^(ntr/2) * V/N ]

    In the functions below:
        temperature:    Temperature (T) in Kelvin
        ntr:            Number of translational degrees of freedom
        mass:           Mass of the molecule in g/mol
        volume:         The volume of the vessel in l
        nmols:          The number of moles in the vessel in mol
        concentration:  Concentration of the system in molar M=mol/l

    The return quanties are in kcal/mol, except for the entropy which is in cal/mol/kelvin, and 
    the logarithm of the partition function, which is a pure number.
*/

/* Compute the translational partition function for ntr degrees of freedom. */
static inline double thermo_tr_lnq(double temperature, int ntr, double mass, double volume, double nmols) {
    double logq_tr  = (0.5*ntr) * log((2.0*PI)*(mass/(AVOGADRO*1000.0))*(BOLTZMANN*temperature)/(PLANCK*PLANCK)) + log((volume/(pow(10.0, ntr)))/(nmols*AVOGADRO));
    return logq_tr;
}

/* Compute the translational (molar) Helmholtz free energy for ntr degrees of freedom. */
static inline double thermo_tr_F(double temperature, int ntr, double mass, double volume, double nmols) {
    return (-J2KCALMOL*BOLTZMANN)*temperature*(thermo_tr_lnq(temperature, ntr, mass, volume, nmols));
}

/* Compute the internal energy associated to ntr translational degrees of freedom. */
static inline double thermo_tr_U(double temperature, int ntr) {
    return (J2KCALMOL*BOLTZMANN)*temperature*ntr*0.5;
}

/* Compute the translational entropy for ntr degrees of freedom. */
static inline double thermo_tr_S(double temperature, int ntr, double mass, double volume, double nmols) {
    return 1000.0 * ( thermo_tr_U(temperature, ntr) - thermo_tr_F(temperature, ntr, mass, volume, nmols) ) / temperature;
}

/* Compute all F, U, S for translation */
static inline void thermo_tr(double temperature, int ntr, double mass, double volume, double nmols, double *LNQ, double *F, double *U, double *S) {
    *LNQ = thermo_tr_lnq(temperature, ntr, mass, volume, nmols);
    *S   = thermo_tr_S(temperature, ntr, mass, volume, nmols);
    *U   = thermo_tr_U(temperature, ntr);
    *F   = thermo_tr_F(temperature, ntr, mass, volume, nmols);
}


/*
    Ideal gas rotations
    -------------------

    Compute thermodynamic quantities for rotational degrees of freedom, for an ideal gas.
    The logarithm of the partition function is:

    log(q_rot) = log ( sqrt(pi)/symmetry ) + 0.5 * sum_i^nrot ( 8 pi^2 I_i kB T / h^2 )

    In the functions below:
        temperature:    Temperature (T) in Kelvin
        nrot:           Number of rotational degrees of freedom
        inertia:        Vector of lenght nrot containing all inertia moments (I_i) in (g/mol)*angstrom^2
        symmetry:       Symmetry number of the molecule

    The return quanties are in kcal/mol, except for the entropy which is in cal/mol/kelvin, and 
    the logarithm of the partition function, which is a pure number.
*/

/* Compute the rotational partition function for one degree of freedom. */
static inline double thermo_rot_lnq_one(double temperature, double inertia) {
    return 0.5*log((((8.0*(PI*PI))*BOLTZMANN)*temperature*(inertia*(1E-23/AVOGADRO)))/(PLANCK*PLANCK));
}

/* Compute the rotational partition function for nrot degrees of freedom. */
static inline double thermo_rot_lnq(double temperature, int nrot, double *inertia, double symmetry) {
    int i;
    double logq_rot = log(sqrt(PI)/symmetry);
    for (i=0; i<nrot; i++) {
        logq_rot += thermo_rot_lnq_one(temperature, inertia[i]);
    }
    return logq_rot;
}

/* Compute the rotational Helmholtz free energy for nrot degrees of freedom. */
static inline double thermo_rot_F(double temperature, int nrot, double *inertia, double symmetry) {
    return (-J2KCALMOL*BOLTZMANN)*temperature*thermo_rot_lnq(temperature, nrot, inertia, symmetry);
}

/* Compute the internal energy associated to nrot rotational degrees of freedom. */
static inline double thermo_rot_U(double temperature, int nrot) {
    return (J2KCALMOL*BOLTZMANN)*temperature*nrot*0.5;
}

/* Compute the rotational entropy for nrot degrees of freedom. */
static inline double thermo_rot_S(double temperature, int nrot, double *inertia, double symmetry) {
    return 1000.0 * ( thermo_rot_U(temperature, nrot) - thermo_rot_F(temperature, nrot, inertia, symmetry) ) / temperature;
}

/* Compute all F, U, S for rotations */
static inline void thermo_rot(double temperature, int nrot, double *inertia, double symmetry, double *LNQ, double *F, double *U, double *S) {
    if (nrot<1.0) {
        *LNQ = 1.0;
        *S   = 0.0;
        *U   = 0.0;
        *F   = 0.0;
    } else {
        *LNQ = thermo_rot_lnq(temperature, nrot, inertia, symmetry);
        *S   = thermo_rot_S(temperature, nrot, inertia, symmetry);
        *U   = thermo_rot_U(temperature, nrot);
        *F   = thermo_rot_F(temperature, nrot, inertia, symmetry);
    }
}


/*
    Ideal gas vibrations (classical)
    --------------------------------

    Compute thermodynamic quantities for vibrational degrees of freedom, for an ideal gas.
    Using classical harmonic oscillator.
    The logarithm of the partition function is:

    log(q_vib,cl) = sum_i^nvib log [ (kB T) / (h v_i) ]

    In the functions below:
        temperature:    Temperature (T) in Kelvin
        nvib:           Number of rotational degrees of freedom
        freq:           Vector of length nvib containing all vibrational frequencies (v_i) in cm-1

    The return quanties are in kcal/mol, except for the entropy which is in cal/mol/kelvin, and 
    the logarithm of the partition function, which is a pure number.
*/

/* Compute the logarithm of the classical vibrational partition function for one vibration */
static inline double thermo_vibcl_lnq_one(double temperature, double freq) {
    return log((BOLTZMANN*temperature)/((PLANCK*LIGHTSPEED*100.0)*freq));
}

/* Compute the logarith of the classical vibration partition function for a set of vibrational modes */
static inline double thermo_vibcl_lnq(double temperature, int nvib, double *freq) {
    int i;
    double logq = 0.0;
    for (i=0; i<nvib; i++) {
        logq += thermo_vibcl_lnq_one(temperature, freq[i]);
    }
    return logq;
}

/* Compute the classical vibration Helmholtz free energy for nvib degrees of freedom. */
static inline double thermo_vibcl_F(double temperature, int nvib, double *freq) {
    return (-J2KCALMOL*BOLTZMANN)*temperature*thermo_vibcl_lnq(temperature, nvib, freq);
}

/* Compute the internal energy associated to nvib classical vibrational degrees of freedom. */
static inline double thermo_vibcl_U(double temperature, int nvib) {
    return (J2KCALMOL*BOLTZMANN)*temperature*nvib;
}

/* Compute the classical vibrational entropy for nvib degrees of freedom. */
static inline double thermo_vibcl_S(double temperature, int nvib, double *freq) {
    return 1000.0 * ( thermo_vibcl_U(temperature, nvib) - thermo_vibcl_F(temperature, nvib, freq) ) / temperature;
}

/* Compute all F, U, S for classical vibrations */
static inline void thermo_vibcl(double temperature, int nvib, double *freq, double *LNQ, double *F, double *U, double *S) {
    *LNQ = thermo_vibcl_lnq(temperature, nvib, freq);
    *S   = thermo_vibcl_S(temperature, nvib, freq);
    *U   = thermo_vibcl_U(temperature, nvib);
    *F   = thermo_vibcl_F(temperature, nvib, freq);
}


/*
    Ideal gas vibrations (quantum)
    --------------------------------

    Compute thermodynamic quantities for vibrational degrees of freedom, for an ideal gas.
    Using quantum harmonic oscillator.
    The logarithm of the partition function is:

    log(q_vib,qm) = sum_i^nvib log [ exp(x)/(1-exp(2x)) ] = - sum_i^nvib log [ 2 sinh(x) ]

    where x = (kB T) / (2 h v_i)

    In the functions below:
        temperature:    Temperature (T) in Kelvin
        nvib :          Number of rotational degrees of freedom
        freq:           Vector of length nvib containing all vibrational frequencies (v_i) in cm-1

    The return quanties are in kcal/mol, except for the entropy which is in cal/mol/kelvin, and 
    the logarithm of the partition function, which is a pure number.
*/

/* Compute the logarithm of the quantum vibrational partition function for one vibration */
static inline double thermo_vibqm_lnq_one(double temperature, double freq) {
    double x = ((PLANCK*LIGHTSPEED*100.0)*freq)/(2.0*BOLTZMANN*temperature);
    return -log(2.0*sinh(x));
}

/* Compute the logarith of the quantum vibration partition function for a set of vibrational modes */
static inline double thermo_vibqm_lnq(double temperature, int nvib, double *freq) {
    int i;
    double logq = 0.0;
    for (i=0; i<nvib; i++) {
        logq += thermo_vibqm_lnq_one(temperature, freq[i]);
    }
    return logq;
}

/* Compute the quantum vibration Helmholtz free energy for nvib degrees of freedom. */
static inline double thermo_vibqm_F(double temperature, int nvib, double *freq) {
    return (-J2KCALMOL*BOLTZMANN)*temperature*thermo_vibqm_lnq(temperature, nvib, freq);
}

/* Compute the internal energy associated to nvib quantum vibrational degrees of freedom. */
static inline double thermo_vibqm_U(double temperature, int nvib, double *freq) {
    int i;
    double x;
    double U=0.0;
    for (i=0; i<nvib; i++) {
        x = ((PLANCK*LIGHTSPEED*100.0)*freq[i])/(2.0*BOLTZMANN*temperature);
        U += x/tanh(x);
    }
    U *= (J2KCALMOL*BOLTZMANN)*temperature;
    return U;
}

/* Compute the quantum vibrational entropy for nvib degrees of freedom. */
static inline double thermo_vibqm_S(double temperature, int nvib, double *freq) {
    return 1000.0 * ( thermo_vibqm_U(temperature, nvib, freq) - thermo_vibqm_F(temperature, nvib, freq) ) / temperature;
}

/* Zero-Point vibrational energy */
static inline double thermo_vibqm_ZPE(int nvib, double *freq) {
    int i;
    double zpe=0.0;
    for (i=0; i<nvib; i++) {
        zpe += (0.5*J2KCALMOL*PLANCK*LIGHTSPEED*100.0)*freq[i];
    }
    return zpe;
}

/* Compute all F, U, S for quantum vibrations */
static inline void thermo_vibqm(double temperature, int nvib, double *freq, double *LNQ, double *F, double *U, double *S, double *ZPE) {
    *LNQ = thermo_vibqm_lnq(temperature, nvib, freq);
    *S   = thermo_vibqm_S(temperature, nvib, freq);
    *U   = thermo_vibqm_U(temperature, nvib, freq);
    *F   = thermo_vibqm_F(temperature, nvib, freq);
    *ZPE = thermo_vibqm_ZPE(nvib, freq);
}


/*
    Solation entropy as by:
        A. J. Garza "Solvation Entropy Made Simple", JCTC 2019.
*/
static inline void solvation_entropy_easysolv(double temperature, double MWs, double Vm, double Vs,
    double density, double acentricity, double permittivity, double expansion,
    double rgyr_m, double rgyr_s, double asa_m, double asa_s,
    double *dStr, double *dSrot, double *Sc_omega, double *Sc_epsilon, double *Sc_alpha) {

    double volume;      /* The corrected volume -- the output! */
    double cbrt_Vm;     /* Cubic root of Vm */
    double cbrt_Vs;     /* Cubic root of Vs */
    double Vfree;       /* Free volume in the solvent (per molecule) */
    double cbrt_Vfree;  /* Cubic root of the free volume */
    double Nc;          /* Average number of accessible cavities */
    double vc;          /* Volume of one solute cavity */
    double cbrt_vc;     /* Cubic root of vc */
    double x;           /* Hopping probability */
    double Nx;          /* Number of sites for hopping per cavity */

    double Ns    = ((AVOGADRO*1000.0*1E-27)*density)/MWs; // Number of molecules per angstrom^3
    double Vone  = 1.0/Ns;  // Angstrom^3 per molecule
    Vfree = Vone - Vs;

    /* Compute vc */
    /*printf("Vs = %lf\n", Vs);
    printf("Vm = %lf\n", Vm);
    printf("Ns = %lf\n", Ns);
    printf("Vone = %lf\n", Vone);
    printf("Vfree = %lf\n", Vfree);*/
    cbrt_Vfree = cbrt(Vfree);
    cbrt_Vm = cbrt(Vm);
    cbrt_vc = cbrt_Vfree + cbrt_Vm;
    vc = cbrt_vc * cbrt_vc * cbrt_vc;
    /*printf("rc = %lf\n", cbrt_vc);
    printf("vc = %lf\n", vc);*/
    double cbrt_vcs = cbrt_Vfree + cbrt(Vs);
    double vcs = cbrt_vcs*cbrt_vcs*cbrt_vcs;
    double rc = cbrt(vc*(3.0/(4.0*PI)));
    double rcs = cbrt(vcs*(3.0/(4.0*PI)));

    /* Compute Nc */
    cbrt_Vs = cbrt(Vs);
    x = fmax( ((cbrt_Vfree*cbrt_Vfree)-(cbrt_Vm*cbrt_Vm)), 0.0 ) / ( (cbrt_Vfree*cbrt_Vfree)+(cbrt_Vs*cbrt_Vs));
    Nx = 4.0 * ( (cbrt_vc*cbrt_vc) / ((cbrt_Vfree*cbrt_Vfree)+(cbrt_Vs*cbrt_Vs)) );
    Nc = 1.0 + Nx * ( x/(1.0-x) );
    volume = Nc*vc*(AVOGADRO*1E-27);
    /*printf("x = %lf\n", x);
    printf("Nx = %lf\n", Nx);
    printf("Nc = %g (~1?)\n", Nc);*/
    //printf("volume = %lf\n", volume);

    /* Translational entropy difference to bring from 1bar gas to liquid confined in volume */
    *dStr = BOLTZMANN*J2KCALMOL*1000.0*log(volume/(BOLTZMANN*AVOGADRO*1000.0*temperature/1E5));
    //printf("dStr = %g\n", dStr);

    // Rotational entropy correction in the liquid
    *dSrot = (3.0*BOLTZMANN*J2KCALMOL*1000.0)*log((rc-rgyr_m)/rc);
    //printf("dSrot = %lf\n", dSrot);

    // Cavity entropy -- omega model
    double Sc_zero = (J2KCALMOL*BOLTZMANN*1000.0) * (Vm/Vs) * log(1.0-Vs*Ns);
    double G = asa_m/asa_s;
    double Sc_omega_dSrot = (3.0*BOLTZMANN*J2KCALMOL*1000.0)*log((rcs-rgyr_s)/rcs);
    *Sc_omega = Sc_zero - G*( (5.365*acentricity*(J2KCALMOL*BOLTZMANN*1000.0)) + Sc_omega_dSrot );
    /*printf("Sc_zero = %lf\n", Sc_zero);
    printf("G = %lf\n", G);
    printf("Sc_omega_dSrot = %lf\n", Sc_omega_dSrot);
    printf("Sc_omega = %lf\n", Sc_omega);*/

    // Cavity entropy -- epsilon model
    double SPT_R = cbrt(Vm/Vs);
    double SPT_y =  ((permittivity-1.0)/(permittivity+2.0))*(3.0/(4.0*PI));
    *Sc_epsilon = -log(1.0-SPT_y) + SPT_R*3.0*(SPT_y/(1.0-SPT_y)) + (SPT_R*SPT_R)*3.0*(SPT_y/(1.0-SPT_y)) + (SPT_R*SPT_R)*(9.0/2.0)*(SPT_y/(1.0-SPT_y))*(SPT_y/(1.0-SPT_y));
    *Sc_epsilon *= -BOLTZMANN*J2KCALMOL*1000.0;
    /*printf("R = %lf\n", SPT_R);
    printf("y = %lf\n", SPT_y);
    printf("Sc_epsilon = %lf\n", Sc_epsilon);*/

    // Cavity entropy -- epsilon-alpha model
    double dfdy = (-(SPT_R*SPT_R)*(6.0*SPT_y+3.0) + 3.0*SPT_R*(SPT_y-1.0) - (SPT_y-1.0)*(SPT_y-1.0)) / ((1.0-SPT_y)*(1.0-SPT_y)*(1.0-SPT_y));
    *Sc_alpha = *Sc_epsilon - ((expansion/1000.0)*temperature)*(BOLTZMANN*J2KCALMOL*1000.0)*SPT_y*dfdy;
    /*printf("dfdy = %lf\n", dfdy);
    printf("Sc_alpha = %lf\n", Sc_alpha);
    printf("Sc_epsilon,alpha = %lf\n", Sc_epsilon + Sc_alpha);*/

}


/*
    Main function to compute everything
*/
double *thermo_compute(double temperature, double energy,
    int ntr, double mass, double volume, double nmols,
    int nrot, double *inertia, double symmetry,
    int nvib, double *freq,
    double solute_volume, double solvent_volume, double solvent_mass, double solvent_density,
    double solvent_acentricity, double solvent_permittivity, double solvent_expansion,
    double rgyr_m, double rgyr_s, double asa_m, double asa_s) {

    /* Vector to store all results */
    double *res;
    res = malloc(THERMO_LAST*sizeof(double));
    if (!res) {
        fprintf(stderr, "Memory allocation failed!\n");
        return NULL;
    }

    /* Energy */
    res[THERMO_LNQ_ELEC] = energy/(J2KCALMOL*BOLTZMANN*temperature);
    res[THERMO_U_ELEC]   = energy;
    res[THERMO_S_ELEC]   = 0.0;
    res[THERMO_F_ELEC]   = energy;

    /* Ideal gas */
    thermo_tr(temperature, ntr, mass, volume, nmols, res+THERMO_LNQ_TR, res+THERMO_F_TR, res+THERMO_U_TR, res+THERMO_S_TR);
    thermo_rot(temperature, nrot, inertia, symmetry, res+THERMO_LNQ_ROT, res+THERMO_F_ROT, res+THERMO_U_ROT, res+THERMO_S_ROT);
    thermo_vibcl(temperature, nvib, freq, res+THERMO_LNQ_VIBCL, res+THERMO_F_VIBCL, res+THERMO_U_VIBCL, res+THERMO_S_VIBCL);
    thermo_vibqm(temperature, nvib, freq, res+THERMO_LNQ_VIBQM, res+THERMO_F_VIBQM, res+THERMO_U_VIBQM, res+THERMO_S_VIBQM, res+THERMO_ZPE);

    /* Sum totals ideal gas */
    res[THERMO_LNQ] = res[THERMO_LNQ_TR] + res[THERMO_LNQ_ROT] + res[THERMO_LNQ_VIBCL] + res[THERMO_LNQ_ELEC];
    res[THERMO_U]   = res[THERMO_U_TR] + res[THERMO_U_ROT] + res[THERMO_U_VIBCL] + res[THERMO_U_ELEC];
    res[THERMO_S]   = res[THERMO_S_TR] + res[THERMO_S_ROT] + res[THERMO_S_VIBCL] + res[THERMO_S_ELEC];
    res[THERMO_F]   = res[THERMO_F_TR] + res[THERMO_F_ROT] + res[THERMO_F_VIBCL] + res[THERMO_F_ELEC];

    /* Solvation entropy based on "Solvation entropy made simple" */
    solvation_entropy_easysolv(temperature, solvent_mass, solute_volume, solvent_volume, solvent_density,
        solvent_acentricity, solvent_permittivity, solvent_expansion, rgyr_m, rgyr_s, asa_m, asa_s,
        res+THERMO_S_EASYSOLV_TR, res+THERMO_S_EASYSOLV_ROT, res+THERMO_S_EASYSOLV_OMEGA, res+THERMO_S_EASYSOLV_EPSILON, res+THERMO_S_EASYSOLV_ALPHA);

    return res;
}


