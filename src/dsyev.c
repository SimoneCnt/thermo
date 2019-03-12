
/*
    Call lapack function DSYEV to diagonalize symmetric matrix.

    Simone Conti 2016-2017
*/

#include <cygtools.h>
#include <thermo.h>

extern void dsyev_(const char* jobz, const char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );

int
mtx_dsyev(int n, double *a, double *w, const char *JOBZ, const char *UPLO) 
{

    #ifndef HAVE_LAPACK
        cyg_logErr("Code compiled without LAPACK support. Impossible to use this functon.");
        return E_FAILURE;
    #else

    int info, lwork;
    double wkopt;
    double* work;

    /* Query and allocate the optimal workspace */
    lwork = -1;
    dsyev_(JOBZ, UPLO, &n, a, &n, w, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work = cyg_malloc(NULL, lwork*cyg_sizeof(double));
    cyg_assert(work!=NULL, E_FAILURE, "Memory allocation failed!");

    /* Solve eigenproblem */
    dsyev_(JOBZ, UPLO, &n, a, &n, w, work, &lwork, &info);

    /* Free workspace */
    free(work);

    /* Check for convergence */
    cyg_assert(info<=0, E_FAILURE, "The algorithm failed to compute eigenvalues.\n");
    return E_SUCCESS;

    #endif
}

