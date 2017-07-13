
/*
    Read a thermo file and save all quantities inside a Thermo structure.

    Copyright (C) 2014-2017 Simone Conti
    Copyright (C) 2015 Université de Strasbourg
*/

#include "cygtools.h"
#include "thermo.h"

/* Read a system from file */
int 
thermo_readthermo(Thermo *A, const char *fname) 
{

    char    *row=NULL, *key, *val;
    int     nr, i;
    FILE    *fp;

    /* Open input config file */
    fp = cyg_fopen(fname, "r");
    cyg_assert(fp!=NULL, E_FAILURE, "Error opening input file.");


    /* Read all elements */
    while (cyg_getline(&row, fp) != -1) {

        /* Skip empty or comment lines */
        if (cyg_isstrempty(row, "# \n\r\t\0")) continue;

        /* Parse all key:val pairs */
        key = strtok(row, "=");
        val = strtok(NULL, "=");
        cyg_assert(val!=NULL, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);

        /* Temperature in kelvin */
        if (strncmp(key, "temperature", 4)==0) {
            nr = sscanf(val, "%lf", &(A->T));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* Number of moles */
        else if (strncmp(key, "nmoles", 4)==0) {
            nr = sscanf(val, "%lf", &(A->n));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* Volume (liters) */
        else if (strncmp(key, "volume", 4)==0) {
            nr = sscanf(val, "%lf", &(A->V));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* Molecular mass in g/mol */
        else if (strncmp(key, "mass", 4)==0) {
            nr = sscanf(val, "%lf", &(A->m));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* Translational degree of freedom */
        else if (strncmp(key, "translations", 4)==0) {
            nr = sscanf(val, "%d", &(A->t));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* Symmetry Number */
        else if (strncmp(key, "sigma", 4)==0) {
            nr = sscanf(val, "%d", &(A->s));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* Rotational degree of freedom and moments of inertia */
        else if (strncmp(key, "rotations", 4)==0) {
            nr = sscanf(val, "%d", &(A->r));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);

            /* Read inertia moments in g/mol/A^2 */
            if (A->r>0) {
                A->I = cyg_malloc(NULL, (A->r)*cyg_sizeof(double));
                cyg_assert(A->I!=NULL, E_FAILURE, "Memory allocation failed!");
                for (i=0; i<A->r; i++) {
                    if (cyg_getline(&row, fp) != -1) {
                        if (cyg_isstrempty(row, "#\n\0")) continue;
                        nr=sscanf(row, "%lf", &(A->I[i]));
                        cyg_assert(nr==1, E_FAILURE, "Impossible to read inertia moment #%d (expected #%d)", i, A->r);
                    }
                }
            }
        }

        /* Vibrational degree of freedom and normal mode frequencies */
        else if (strncmp(key, "vibrations", 4)==0) {
            nr = sscanf(val, "%d", &(A->v));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);

            /* Read vibrational modes in cm-1 */
            A->nu = cyg_malloc(NULL, (A->v)*cyg_sizeof(double));
            cyg_assert(A->nu!=NULL, E_FAILURE, "Memory allocation failed!");
            for (i=0; i<A->v; i++) {
                if (cyg_getline(&row, fp) != -1) {
                    if (cyg_isstrempty(row, "#\n\0")) continue;
                    nr=sscanf(row, "%lf", &(A->nu[i]));
                    cyg_assert(nr==1, E_FAILURE, "Impossible to read vibration #%d (expected #%d)", i, A->v);
                }
            }
        }

        /* Energy in kcal/mol */
        else if (strncmp(key, "energy", 4)==0) {
            nr = sscanf(val, "%lf", &(A->E));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* Hessian matrix file */
        else if (strncmp(key, "hessian", 4)==0) {
            char tmpstr[128];
            nr = sscanf(val, "%127s", tmpstr);
            A->hessfile = strdup(tmpstr);
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* Unknown Keyword */
        else {
            cyg_logErr("Unknown keyword <%s>", key);
            return E_FAILURE;
        }
    }

    /* Clean and return */
    fclose(fp);
    free(row);
    return E_SUCCESS;
}

