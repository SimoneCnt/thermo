
/*
    Read a thermo file and save all quantities inside a Thermo structure.

    Copyright (C) 2014-2019 Simone Conti
    Copyright (C) 2015 Université de Strasbourg
*/

#include <cygtools.h>
#include <thermo.h>

/* Read a system from file */
int 
thermo_readthermo(Thermo *A, const char *fname) 
{

    char    *row=NULL, *key, *val;
    char    unit[8];
    int     nr, i;
    double  tmpd, concentration;
    FILE    *fp;
    bool    convert_density_to_volume=false;

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

        /* Concentration */
        else if (strncmp(key, "concentration", 13)==0) {
            memset(unit, '\0', sizeof(unit));
            nr = sscanf(val, "%lf", &concentration);
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
            if (strstr(val, "unit")!=NULL) {
                nr = sscanf(val, " %*f %*s %7s", unit);
                cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s> while reading unit", val, key);
                if (strncmp(unit, "M", 1)==0) {
                    fprintf(fpout, "Found unit <%s> for concentration\n", unit);
                    A->V = concentration;
                    A->n = 1.0;
                } else if (strncmp(unit, "g/ml", 4)==0) {
                    fprintf(fpout, "Found unit <%s> for concentration\n", unit);
                    A->V = 1.0/(1000.0*concentration);
                    convert_density_to_volume = true;
                } else {
                    cyg_logErr("Impossible to understand unit <%s> for concentration. Possible values are M (default) or g/ml\n", unit);
                    return E_FAILURE;
                }
            }
        }

        /* Pressure (atm) */
        else if (strncmp(key, "pressure", 4)==0) {
            nr = sscanf(val, "%lf", &(A->pressure));
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
            memset(unit, '\0', sizeof(unit));
            nr = sscanf(val, "%d", &(A->r));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
            if (strstr(val, "unit")!=NULL) {
                nr = sscanf(val, " %*d %*s %7s", unit);
                cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s> while reading unit", val, key);
                if (strncmp(unit, "K", 1)==0 || strncmp(unit, "gmolA2", 6)==0 || strncmp(unit, "cm-1", 4)==0) {
                    fprintf(fpout, "Found unit <%s> for rotations\n", unit);
                } else {
                    cyg_logErr("Impossible to understand unit <%s> for rotations. Possible values are gmolA2 (default) or K\n", unit);
                    return E_FAILURE;
                }
            }
            /* Read inertia moments */
            if (A->r>0) {
                A->I = cyg_malloc(NULL, (A->r)*cyg_sizeof(double));
                cyg_assert(A->I!=NULL, E_FAILURE, "Memory allocation failed!");
                for (i=0; i<A->r; i++) {
                    if (cyg_getline(&row, fp) != -1) {
                        if (cyg_isstrempty(row, "#\n\0")) continue;
                        nr=sscanf(row, "%lf", &tmpd);
                        cyg_assert(nr==1, E_FAILURE, "Impossible to read inertia moment #%d (expected #%d)", i, A->r);
                        if (strncmp(unit, "K", 1)==0) {
                            A->I[i] = thermo_kelvin2inertia(tmpd);
                        } else if (strncmp(unit, "cm-1", 4)==0) {
                            A->I[i] = thermo_freq2inertia(tmpd);
                        } else {
                            A->I[i] = tmpd;
                        }
                    }
                }
            }
        }

        /* Vibrational degree of freedom and normal mode frequencies */
        else if (strncmp(key, "vibrations", 4)==0) {
            memset(unit, '\0', sizeof(unit));
            nr = sscanf(val, "%d", &(A->v));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
            if (strstr(val, "unit")!=NULL) {
                nr = sscanf(val, " %*d %*s %7s", unit);
                cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s> while reading unit", val, key);
                if (strncmp(unit, "K", 1)==0 || strncmp(unit, "cm-1", 6)==0) {
                    fprintf(fpout, "Found unit <%s> for vibrations\n", unit);
                } else {
                    cyg_logErr("Impossible to understand unit <%s> for vibrations. Possible values are cm-1 (default) or K\n", unit);
                    return E_FAILURE;
                }
            }
            /* Read vibrational modes */
            if (A->v<1) continue;
            A->nu = cyg_malloc(NULL, (A->v)*cyg_sizeof(double));
            cyg_assert(A->nu!=NULL, E_FAILURE, "Memory allocation failed!");
            for (i=0; i<A->v; i++) {
                if (cyg_getline(&row, fp) != -1) {
                    if (cyg_isstrempty(row, "#\n\0")) continue;
                    nr=sscanf(row, "%lf", &tmpd);
                    cyg_assert(nr==1, E_FAILURE, "Impossible to read vibration #%d (expected #%d)", i, A->v);
                    if (strncmp(unit, "K", 1)==0) {
                        A->nu[i] = thermo_kelvin2cm(tmpd);
                    } else {
                        A->nu[i] = tmpd;
                    }
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

        /* van der Waals volume of the solute [A^3] */
        else if (strncmp(key, "solute_volume", 13)==0) {
            nr = sscanf(val, "%lf", &(A->solute_volume));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* van der Waals volume of the solute [A^3] */
        else if (strncmp(key, "solvent_volume", 14)==0) {
            nr = sscanf(val, "%lf", &(A->solvent_volume));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* solvent molecular weight [g/mol] */
        else if (strncmp(key, "solvent_mass", 12)==0) {
            nr = sscanf(val, "%lf", &(A->solvent_mass));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* density of the solvent [g/ml] */
        else if (strncmp(key, "solvent_density", 15)==0) {
            nr = sscanf(val, "%lf", &(A->density));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* acentricity of the solvent */
        else if (strncmp(key, "solvent_acentricity", 19)==0) {
            nr = sscanf(val, "%lf", &(A->acentricity));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* permittivity of the solvent */
        else if (strncmp(key, "solvent_permittivity", 20)==0) {
            nr = sscanf(val, "%lf", &(A->permittivity));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* isobaric thermcal expansion coefficient of the solvent [10^-3/K] */
        else if (strncmp(key, "solvent_expansion", 17)==0) {
            nr = sscanf(val, "%lf", &(A->thermal_expansion));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* radius of gyration of the solute [A] */
        else if (strncmp(key, "rgyr", 4)==0) {
            nr = sscanf(val, "%lf", &(A->rgyr_m));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* radius of gyration of the solvent [A] */
        else if (strncmp(key, "solvent_rgyr", 12)==0) {
            nr = sscanf(val, "%lf", &(A->rgyr_s));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* sasa of the solute [A^2] */
        else if (strncmp(key, "sasa", 4)==0) {
            nr = sscanf(val, "%lf", &(A->asa_m));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* sasa of the solvent [A^2] */
        else if (strncmp(key, "solvent_sasa", 12)==0) {
            nr = sscanf(val, "%lf", &(A->asa_s));
            cyg_assert(nr==1, E_FAILURE, "Invalid value <%s> for key <%s>", val, key);
        }

        /* Unknown Keyword */
        else {
            cyg_logErr("Unknown keyword <%s>", key);
            return E_FAILURE;
        }
    }

    /* Validate nmols, volume, concentration, and pressure */
    if (A->pressure>0.0) {
        A->n = 1.0;
        A->V = (1000.0*CNS_kB*CNS_NA) * A->T / (A->pressure*101325.0);
    }
    if (convert_density_to_volume) {
        A->V *= A->m;
    }

    /* Clean and return */
    fclose(fp);
    free(row);
    return E_SUCCESS;
}

