/******************************************************************************
 *
 *  thermo/readthermo.c
 *  Read a thermo file and save all quantities inside a Thermo structure
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
 * @brief Read a thermo file and save all quantities inside a Thermo structure
 *
 * Read a thermo file and save all quantities inside a Thermo structure
 * 
 * @param[in,out] A Pointer to an initialized @c Thermo structure
 * @param[in] fname Input thermo filename
 * @retval EXIT_SUCCESS if everithing is ok, EXIT_FAILURE otherwise
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermo.h"

/* Read a system from file */
int 
thermo_readthermo(Thermo *A, const char *fname) 
{

    char    *tmp, row[180], c[2];
    int     nr, i;
    double  val;
    FILE    *fp;

    /* Open input config file */
    fp = fopen(fname, "r");
    if (fp==NULL) {
            fprintf(stderr,"\n   ERROR!! Is not possible to open %s in mode r. Quit.\n\n", fname);
            exit(EXIT_FAILURE);
    }

    /* Read all elements */
    while ((tmp=fgets(row, 180, fp))!=NULL) {
        nr=sscanf(row,"%1s %lf", c, &val);
        if (nr<=0) {            /* Skip empty lines */
            continue;
        } else if (c[0]=='#') { /* Skip comments */
            continue;
        } else if (nr==1) {
            fprintf(stderr,"\n   ERROR!! I cannot understand your input file at \n");
            fprintf(stderr,"      %s", row);
            fprintf(stderr,"   Quitting.\n\n");
            exit(EXIT_FAILURE);
        }
        switch (c[0]) {
            case 'T':           /* Temperature in kelvin */
                A->T = val;
                break;
            case 'n':           /* Number of moles */
                A->n = val;
                break;
            case 'V':           /* Volume (liters) */
                A->V = val;
                break;
            case 'm':           /* Molecular mass in g/mol */
                A->m = val; 
                break;
            case 't':           /* Translational degree of freedom */
                A->t = val;
                break;
            case 's':           /* Symmetry Number */
                A->s = val;
                break;
            case 'r':           /* Rotational degree of freedom and moments of inertia */
                A->r = val;
                A->I=(double*)malloc(A->r*sizeof(double));
                for (i=0; i<A->r; i++) {    /* all inertia moments in g/mol/A^2 */
                    if ((tmp=fgets(row, 80, fp))!=NULL) {
                        if (row[0]!='#' && row[0]!='\n') {	
                            if ((nr=sscanf(row,"%lf",&val))!=1) {
                                fprintf(stderr,"\n   ERROR r !! fscanf return %d instead of 1\n",nr);
                                exit(EXIT_FAILURE);
                            }
                            A->I[i] = val;
                        }
                    }
                }
                break;
            case 'v':           /* Vibrational degree of freedom and normal mode frequencies */
                A->v = val;
                A->nu=(double*)malloc(A->v*sizeof(double));
                for (i=0; i<A->v; i++) {    /* all vibrational modes in cm-1 */
                    if ((tmp=fgets(row, 80, fp))!=NULL) {
                        if (row[0]!='#' && row[0]!='\n') {	
                            if ((nr=sscanf(row,"%lf",&val))!=1) {
                                fprintf(stderr,"\n  ERROR v !! fscanf return %d instead of 1\n",nr);
                                exit(EXIT_FAILURE);
                            }
                            A->nu[i]=val;
                        }
                    }
                }
                break;		
            case 'E':           /* Energy in kcal/mol */
                A->E = val;
                break;
            default:
                printf("Warning: unknown option %s\n",c);
                break;
        }
    }
    fclose(fp);
    return 0;
}

