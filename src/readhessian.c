
/*
    Read hessian matrix printed by CHARMM.

    Simone Conti 2016-2017
*/

#include <cygtools.h>
#include <thermo.h>

int 
thermo_readhessian(Thermo *A)
{

    fprintf(fpout, "Reading hessian file <%s>...\n", A->hessfile);

    char *filename=A->hessfile;
    FILE *fp=cyg_fopen(filename, "r");
    cyg_assert(fp!=NULL, E_FAILURE, "Impossible to open <%s>!\n", filename);

    int i, j, nf, nat, nat3;
    double tmp;

    /* Read the number of atoms */
    nf=fscanf(fp, "%d", &nat);
    cyg_assert(nf==1, E_FAILURE, "No number of atoms read!\n");
    cyg_assert(nat>1, E_FAILURE, "Number of atoms = %d (<2)\n", nat);
    nat3 = nat*3;

    /* Read (and discard) the energy */
    nf=fscanf(fp, "%lf", &tmp);
    cyg_assert(nf==1, E_FAILURE, "Unexpected end of file.\n");

    /* Read (and discard) the gradient */
    for (i=0; i<nat3; i++) {
        nf=fscanf(fp, "%lf", &tmp);
        cyg_assert(nf==1, E_FAILURE, "Unexpected end of file.\n");
    }

    /* Read the hessian */
    A->hessian = cyg_malloc(A->hessian, nat3*nat3*cyg_sizeof(double));
    cyg_assert(A->hessian!=NULL, E_FAILURE, "Memory allocation failed!");
    for (i=0; i<nat3; i++) {
        for (j=i; j<nat3; j++) {
            nf=fscanf(fp, "%lf", &tmp);
            cyg_assert(nf==1, E_FAILURE, "Unexpected end of file.\n");
            A->hessian[i+j*nat3] = A->hessian[j+i*nat3] = tmp;
        }
    }
    A->natoms = nat;

    fclose(fp);
    return E_SUCCESS;

}

