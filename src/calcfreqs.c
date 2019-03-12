
/*
    Diagonalize a hessian matrix and calculate vibrational frequences.

    Simone Conti 2016-2017
*/

#include <cygtools.h>
#include <thermo.h>

int
thermo_calcfreqs(Thermo *A)
{

    fprintf(fpout, "Diagonalizing hessian matrix and calculating frequencies...\n");

    int i, nat3, skip;
    double *eival;

    nat3 = A->natoms*3;

    /* Eigenvalues */
    eival = cyg_malloc(NULL, nat3*cyg_sizeof(double));
    cyg_assert(eival!=NULL, E_FAILURE, "Memory allocation failed!");

    /* Frequencies */
    A->nu = cyg_malloc(A->nu, nat3*cyg_sizeof(double));
    cyg_assert(A->nu!=NULL, E_FAILURE, "Memory allocation failed!");

    /* Diagonalize matrix */
    mtx_dsyev(nat3, A->hessian, eival, "V", "L");

    /* Convert eigenvalues to frequencies */
    skip = A->t + A->r;
    fprintf(fpout, "Number of atoms: %d\n", A->natoms);
    fprintf(fpout, "Total number of degrees of freedom: %d\n", nat3);
    fprintf(fpout, "Skipping %d for translations and rotations.\n", skip);
    fprintf(fpout, "Obtained %d vibrational modes.\n", nat3-skip);
    double cvtfrq=2045.5/(2.99793*6.28319);
    for (i=0; i<nat3-skip; i++) {
        A->nu[i] = cvtfrq*sqrt(fabs(eival[i+skip]))*copysign(1.0, eival[i+skip]);
    }
    A->v = nat3-skip;

    /* Clean and return */
    free(eival);
    return E_SUCCESS;
}

