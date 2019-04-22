
/*
    Thermo

    Copyright (C) 2014-2019 Simone Conti
    Copyright (C) 2015-2016 Université de Strasbourg
*/

#include <cygtools.h>
#include <thermo.h>


/* Functions defined at the end of this file */
static void version(void);  /* Print version info */
static void version2(void);  /* Print version info -- more system specifics*/
static void usage(void);    /* Print the usage of the software */
static void help(void);     /* Print some help */

FILE *fpout=NULL;

/* Main */
int 
main(int argc, char *argv[]) 
{

    /* Declare used variables */
    int hasA=0, hasB=0, hasStechio=0, nA, nB, nr, cumul=0, vdos=0, ret;
    char *nameA=NULL, *nameB=NULL;
    char *outfile=NULL;
    bool raw_output = false;
    fpout = stderr;

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
        {"out",     required_argument, 0, 'o'},
        {"raw",     no_argument,       0, 'r'},
        {"stechio", required_argument, 0, 's'},
        {"cumul",   no_argument,       0, 'c'},
        {"vdos",    no_argument,       0, 'd'},
        {"dnu",     required_argument, 0, 'n'},
        {"version", no_argument,       0, 'v'},
        {"help",    no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };


    /* Parse command line options */
    while (1) {
        c = getopt_long_only(argc, argv, "A:B:o:rs:cdn:vh", long_options, &option_index);

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

            case 'o': /* Output file */
                outfile = optarg;
                break;

            case 'r': /* raw output format */
                raw_output = true;
                break;

            case 's': /* Stechiometric coefficients */
                nr = sscanf(optarg, "%d:%d", &nA, &nB);
                if (nr!=2) {
                    version();
                    fprintf(stderr, "Error parsing --stechio option! I cannot find two coefficients in a:b form!\n\n");
                    usage();
                    return EXIT_FAILURE;
                } else if (nA<1 || nB<1) {
                    version();
                    fprintf(stderr, "Error parsing --stechio option! One (or both) of the two coefficients are less than one!\n\n");
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
                break;

            case 'n': /* Accuracy in vibration hystograms */
                sscanf(optarg, "%lf", &A.dnu);
                A.nu_np = (int)lrint(ceil(4000.0/A.dnu));
                B.dnu   = A.dnu;
                B.nu_np = A.nu_np;
                D.dnu   = A.dnu;
                D.nu_np = A.nu_np;
                break;

            case 'v': /* Version */
                version();
                version2();
                return EXIT_SUCCESS;
                break;

            case 'h': /* Help */
                version();
                help();
                return EXIT_SUCCESS;
                break;

            case '?': /* Error in command line */
                version();
                fprintf(stderr, "Error parsing command line options. Exiting.\n\n");
                usage();
                return EXIT_FAILURE;
                break;

            default: /* Unknown error */
                version();
                fprintf(stderr, "Unexpected error in command line parsing. Exiting.\n\n");
                usage();
                return EXIT_FAILURE;
        }
    }

    /* Extraneous options */
    if (optind < argc) {
        version();
        fprintf(stderr, "Error! Extraneous options found in command line parsing: ");
        while (optind < argc) {
            fprintf(stderr, "%s ", argv[optind++]);
        }
        fprintf(stderr, "\n\n");
        usage();
        return EXIT_FAILURE;
    }

    /* Open outfile for writing */
    if (outfile) {
        fpout = fopen(outfile, "w");
        if (!fpout) {
            version();
            fprintf(stderr, "ERROR! Impossible to open file <%s> for writing!\n", outfile);
            return EXIT_FAILURE;
        }
    } else {
        fpout = stdout;
    }

    /* Print version */
    if (!raw_output) version();

    /* Check if you gave at least A or B */
    if (!hasA && !hasB) {
        fprintf(stderr, "Error! You did not specified neither A nor B!\n");
        return EXIT_FAILURE;
    }

    /* Command line parsing went ok. Can continue. */

    /* Work with mol A */
    if (hasA) {
        char *sname = strrchr(nameA, '/'); if (sname==NULL) sname=nameA; else sname++;
        if (!(raw_output && !hasStechio)) fprintf(fpout, "\nMolecule A: <%s>\
                \n---------------------------------------------\n\n", sname);
        ret = thermo_readthermo(&A, nameA);
        cyg_assert(ret==E_SUCCESS, E_FAILURE, "Failing reading thermo input file <%s>", nameA);
        if (A.hessfile) {
            thermo_readhessian(&A);
            thermo_calcfreqs(&A);
        }
        thermo_printconfig(&A, raw_output);
        thermo_calcthermo(&A);
        thermo_printthermo(&A,0, raw_output);
        if (cumul) thermo_cumulvib(&A, "cumul_A");
        if (vdos)  thermo_vdos(&A, "vdos_A.dat");
    }

    /* Work with mol B */
    if (hasB) {
        char *sname = strrchr(nameB, '/'); if (sname==NULL) sname=nameB; else sname++;
        if (!(raw_output && !hasStechio)) fprintf(fpout, "\nMolecule B: <%s>\
                \n---------------------------------------------\n\n", sname);
        ret = thermo_readthermo(&B, nameB);
        cyg_assert(ret==E_SUCCESS, E_FAILURE, "Failing reading thermo input file <%s>", nameB);
        if (B.hessfile) {
            thermo_readhessian(&B);
            thermo_calcfreqs(&B);
        }
        thermo_printconfig(&B, raw_output);
        thermo_calcthermo(&B);
        thermo_printthermo(&B,0, raw_output);
        if (cumul) thermo_cumulvib(&B, "cumul_B");
        if (vdos)  thermo_vdos(&B, "vdos_B.dat");
    }

    /* Evaluate difference in reaction */
    if (hasA && hasB && hasStechio) {
        fprintf(fpout, "\nDifferences for the reaction %dA <-> %dB\
                \n---------------------------------------------\n\n", nA, nB);
        thermo_diffthermo(&A, &B, nA, nB, &D);
        thermo_printthermo(&D,1, raw_output);
        if (cumul) thermo_cumulvib(&D, "cumul_D");
    }

    /* Cleaning */
    thermo_delete(&A);
    thermo_delete(&B);
    thermo_delete(&D);

    return 0;
}

void version() {
    fprintf(fpout, "\n");
    fprintf(fpout, "    Thermo 2.0\n");
    fprintf(fpout, "    ==========\n");
    fprintf(fpout, "\n");
    fprintf(fpout, "Copyright (C) 2014-2017-2019 Simone Conti\n");
    fprintf(fpout, "Copyright (C) 2015-2016 Université de Strasbourg\n");
    fprintf(fpout, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\n");
    fprintf(fpout, "This is free software: you are free to change and redistribute it.\n");
    fprintf(fpout, "There is NO WARRANTY, to the extent permitted by law.\n");
    fprintf(fpout, "\n");
    fprintf(fpout, "Written by Simone Conti.\n");
    fprintf(fpout, "\n");
}

void version2() {
    fprintf(fpout, "GIT version: %s\n", GIT_VERSION);
    fprintf(fpout, "Compiled on %s using the %s C compiler v%s \n", BUILD_DATE, CC_ID, CC_VERSION);
    fprintf(fpout, "    on a %s (%s) machine\n", SYSTEM_GEN, SYSTEM_PROC);
    fprintf(fpout, "\n");
}

void usage() {
    fprintf(fpout, "Usage: thermo [OPTION]...\n");
    fprintf(fpout, "\n");
    fprintf(fpout, "Options:\n");
    fprintf(fpout, "   -A, --A        fname   Input thermo file for the molecule A\n");
    fprintf(fpout, "   -B, --B        fname   Input thermo file for the molecule B\n");
    fprintf(fpout, "   -o, --out      fname   Output file\n");
    fprintf(fpout, "   -s, --stechio  a:b     Stechiometric coefficients for the reaction aA<->bB\n");
    fprintf(fpout, "   -c, --cumul    fname   Print the cumulative vibrational chemical potential\n");
    fprintf(fpout, "   -d, --vdos     fname   Print the vibrational density of state\n");
    fprintf(fpout, "   -n, --dnu      real    Accuracy in the calculation of the vibration hystograms\n");
    fprintf(fpout, "   -h, --help             Show this help and exit\n");
    fprintf(fpout, "   -v, --version          Print version information and exit\n");
    fprintf(fpout, "\n");
}

void help() {
    fprintf(fpout, "thermo is a wonderfull software :D\n");
    fprintf(fpout, "\n");
    usage();
    fprintf(fpout, "\n");
    fprintf(fpout, "Examples:\n");
    fprintf(fpout, "  See examples/ directory for working examples.\n");
    fprintf(fpout, "\n");
    fprintf(fpout, "Report bugs to <https://github.com/SimoneCnt/thermo/issues> \n");
    fprintf(fpout, "  or directly to <simonecnt@gmail.com>.\n");
    fprintf(fpout, "\n");
}

