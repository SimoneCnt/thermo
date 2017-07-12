
/*
    Some macro and functions...

    Simone Conti 2017
*/

#ifndef _CYGNE_TOOLS_H_
#define _CYGNE_TOOLS_H_

/* Include "standard" libraries */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <sys/time.h>
#include <getopt.h>
#include <errno.h>
#include <string.h>

#ifdef HAVE_QUADMATH
    #include <quadmath.h>
    typedef  __float128 quad;
#endif

#ifdef HAVE_THREADS
    #include <pthread.h>
#endif

#ifndef M_PI
#define M_PI        3.14159265358979323846  /* pi */
#endif

/* Defined by complex.h */
#undef I

/* Return codes */
#define E_SUCCESS EXIT_SUCCESS
#define E_FAILURE EXIT_FAILURE

/* Define CYG_FILE, CYG_FUNC and CYG_LINE to use in error messages. */
#if __STDC_VERSION__ < 199901L
# if __GNUC__ >= 2
#  define __func__ __FUNCTION__
# else
#  define __func__ "<unknown>"
# endif
#endif
#define CYG_FILE __FILE__
#define CYG_FUNC __func__
#define CYG_LINE __LINE__

/* Print current function, position and line. */
#define cyg_printsource(fp) fprintf(fp, "%s (file %s, line %d):\n", CYG_FUNC, CYG_FILE, CYG_LINE);

/* Print an error. */
#define cyg_logErr(...) {\
    fprintf(stderr, "\n"); \
    cyg_printsource(stderr); \
    fprintf(stderr, "    ERROR!! "); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
}

/* Check an expression. If false print message and return error code. */
#define cyg_assert(expr, ret, ...) {\
    if (!(expr)) {\
        cyg_logErr("Assert failed! %s", #expr);\
        fprintf(stderr, __VA_ARGS__);\
        fprintf(stderr, "\n");\
        return ret;\
    }\
}

/* Sizeof that takes int instead of size_t... */
#define cyg_sizeof(x) ((int)sizeof(x))

/* Open a file */
static inline FILE *cyg_fopen(const char *filename, const char *mode) {
    cyg_assert(filename!=NULL, NULL, "You didn't specify any filename");
    cyg_assert(mode!=NULL, NULL, "You didn't specify an opening modality");
    FILE *fp = fopen(filename, mode);
    cyg_assert(fp!=NULL, NULL, "Error opening %s in mode %s: %s", filename, mode, strerror(errno));
    return fp;
}

/* (re)Allocate/free memory. *p must be null or a valid pointer from cyg_malloc */
static inline void *cyg_malloc(void *p, int n) {
    if (n <= 0) { free(p); return NULL; }
    void *pp = realloc(p, (size_t)n);
    cyg_assert(pp!=NULL, NULL, "Memory allocation failed fo %d bytes! Not enough memory!", n);
    return pp;
}

/* Return true if str[0] is equal to any char in vals */
static inline bool cyg_isstrempty(char *str, const char *vals) {
    if (vals==NULL) return false;
    size_t i;
    for (i=0; i<strlen(vals); i++) {
        if (str[0]==vals[i]) return true;
    }
    return false;
}    

/* Read a line from a stream */
int cyg_getline (char **lineptr, FILE *stream);

#endif

