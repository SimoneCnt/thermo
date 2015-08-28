/******************************************************************************
 *
 *  thermo/thermo.c
 *
 *  Copyright (C) 2014, 2015 Simone Conti
 *  Copyright (C) 2015 Université de Strasbourg
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include "thermo.h"

/* Functions defined at the end of this file */
static void version(void);  /* Print version info */
static void usage(void);    /* Print the usage of the software */
static void help(void);     /* Print some help */

/* Main */
int 
main(int argc, char *argv[]) 
{

    /* Declare used variables */
    int hasA=0, hasB=0, hasStechio=0, nA, nB, nr, cumul=0, vdos=0;
    char *nameA, *nameB, *namevdos;

    /* Define and initialize Thermo structures */
    Thermo A, B, D;
    thermo_init(&A);
    thermo_init(&B);
    thermo_init(&D);

    /* Getopt variables */
    int option_index=0, c;
    static struct option long_options[] = {
        {"A",       required_argument, 0, 'A'},
        {"B",       required_argument, 0, 'B'},
        {"stechio", required_argument, 0, 's'},
        {"cumul",   required_argument, 0, 'c'},
        {"vdos",    required_argument, 0, 'd'},
        {"dnu",     required_argument, 0, 'n'},
        {"version", no_argument,       0, 'v'},
        {"help",    no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    version();

    /* Parse command line options */
    while (1) {
        c = getopt_long_only(argc, argv, "A:B:s:cd:n:vh", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1) break;

        switch (c) {

            case 'A': /* Molecule A */
                hasA  = 1;
                nameA = optarg;
                break;

            case 'B': /* Molecule B */
                hasB  = 1;
                nameB = optarg;
                break;

            case 's': /* Stechiometric coefficients */
                nr = sscanf(optarg, "%d:%d", &nA, &nB);
                if (nr!=2) {
                    printf("Error parsing --stechio option! I cannot find two coefficients in a:b form!\n\n");
                    usage();
                    return EXIT_FAILURE;
                } else if (nA<1 || nB<1) {
                    printf("Error parsing --stechio option! One (or both) of the two coefficients are less than one!\n\n");
                    usage();
                    return EXIT_FAILURE;
                } else {
                    hasStechio = 1;
                }
                break;

            case 'c': /* Cumulative vibrational chemical potential */
                cumul=1;
                break;

            case 'd': /* Vibrational density of states */
                vdos=1;
                namevdos = optarg;
                break;

            case 'n': /* Accuracy in vibration hystograms */
                sscanf(optarg, "%lf", &A.dnu);
                A.nu_np = lrint(ceil(4000.0/A.dnu));
                B.dnu   = A.dnu;
                B.nu_np = A.nu_np;
                D.dnu   = A.dnu;
                D.nu_np = A.nu_np;
                break;

            case 'v': /* Version */
                return EXIT_SUCCESS;
                break;

            case 'h': /* Help */
                help();
                return EXIT_SUCCESS;
                break;

            case '?': /* Error in command line */
                fprintf(stderr, "Error parsing command line options. Exiting.\n\n");
                usage();
                return EXIT_FAILURE;
                break;

            default: /* Unknown error */
                printf("Unexpected error in command line parsing. Exiting.\n\n");
                usage();
                return EXIT_FAILURE;
        }
    }

    /* Extraneous options */
    if (optind < argc) {
        printf ("Error! Extraneous options found in command line parsing: ");
        while (optind < argc) {
            printf("%s ", argv[optind++]);
        }
        printf("\n\n");
        usage();
        return EXIT_FAILURE;
    }

    /* Check if you gave at least A or B */
    if (!hasA && !hasB) {
        printf("Error! You did not specified neither A nor B!\n");
        return EXIT_FAILURE;
    }

    /* Command line parsing went ok. Can continue. */

    /* Work with mol A */
    if (hasA) {
        printf("\nMolecule A: <%s>\
                \n---------------------------------------------\n\n", nameA);
        thermo_readthermo(&A, nameA);
        thermo_printconfig(&A);
        thermo_calcthermo(&A);
        thermo_printthermo(&A,0);
        if (cumul) thermo_cumulvib(&A, "fvib");
        if (vdos)  {
            thermo_vdos(&A, namevdos);
            thermo_vdosfvib(&A,"fvib.vdos.dat");
        }
    }

    /* Work with mol B */
    if (hasB) {
        printf("\nMolecule B: <%s>\
                \n---------------------------------------------\n\n", nameB);
        thermo_readthermo(&B, nameB);
        thermo_printconfig(&B);
        thermo_calcthermo(&B);
        thermo_printthermo(&B,0);
        if (cumul) thermo_cumulvib(&B, "molB.cumul");
        if (vdos)  thermo_vdos(&B, namevdos);
    }

    /* Evaluate difference in reaction */
    if (hasA && hasB && hasStechio) {
        printf("\nDifferences for the reaction %dA <-> %dB\
                \n---------------------------------------------\n\n", nA, nB);
        thermo_diffthermo(&A, &B, nA, nB, &D);
        thermo_printthermo(&D,1);
    }

    /* Cleaning */
    thermo_delete(&A);
    thermo_delete(&B);
    thermo_delete(&D);

    printf("\n");
    return 0;
}



/* Initial concentration scan */
/*int main(int argc, char *argv[]) {
    thermo A;
    cyg_thermo_init(&A);
    cyg_thermo_read(&A, argv[1]);
    double C0, mu;
    for (C0=1e-12; C0<10; C0*=1.2) {
        A.C = C0;
        mu = print_all(&A);
        printf(" %+e %lf \n" , C0, mu);
    }
    free(A.nu);
    free(A.I);
    return 0;
}*/


void version() {
    printf("\n");
    printf("    Thermo 1.0\n");
    printf("    ==========\n");
    printf("\n");
    printf("Copyright (C) 2014, 2015 Simone Conti\n");
    printf("Copyright (C) 2015 Université de Strasbourg\n");
    printf("License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\n");
    printf("This is free software: you are free to change and redistribute it.\n");
    printf("There is NO WARRANTY, to the extent permitted by law.\n");
    printf("\n");
    printf("Written by Simone Conti.\n");
    printf("\n");
}

void usage() {
    printf("Usage: thermo [OPTION]...\n");
    printf("\n");
    printf("Options:\n");
    printf("   -A, --A        fname   Input thermo file for the molecule A\n");
    printf("   -o, --B        fname   Input thermo file for the molecule B\n");
    printf("   -s, --stechio  a:b     Stechiometric coefficients for the reaction aA<->bB\n");
    printf("   -c, --cumul    fname   Print the cumulative vibrational chemical potential\n");
    printf("   -d, --vdos     fname   Print the vibrational density of state\n");
    printf("   -n, --dnu      real    Accuracy in the calculation of the vibration hystograms\n");
    printf("   -h, --help             Show this help and exit\n");
    printf("   -v, --version          Print version information and exit\n");
    printf("\n");
}

void help() {
    printf("thermo is a wonderfull software :D\n");
    printf("\n");
    usage();
    printf("\n");
    printf("Examples:\n");
    printf("  See test/ directory for working examples.\n");
    printf("\n");
    printf("Report bugs to <https://github.com/SimoneCnt/thermo/issues> \n");
    printf("  or directly to <simonecnt@gmail.com>.\n");
    printf("\n");
}

