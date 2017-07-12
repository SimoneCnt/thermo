 
/*
    Read up to (and including) a new line terminator from an input stream.

    Copyright (C) 1993 Free Software Foundation, Inc.
    Copyright (C) 2014-2017 Simone Conti
*/

#include <assert.h>
#include "cygtools.h"

int cyg_getline (char **lineptr, FILE *stream) {

    int nchars_avail;		/* Allocated but unused chars in *LINEPTR.  */
    char *read_pos;		/* Where we're reading into *LINEPTR. */
    int ret;
    int nn=1;
    int *n=&nn;

    /* Always add at least this many bytes when extending the buffer.  */
    int min_chunk = 64;

    if (!lineptr || !n || !stream) {
        errno = EINVAL;
        return -1;
    }

    if (!*lineptr) {
        *n = min_chunk;
        *lineptr = malloc ((size_t)*n);
        if (!*lineptr) {
            errno = ENOMEM;
            return -1;
        }
    }

    nchars_avail = *n;
    read_pos = *lineptr;

    for (;;) {
        int save_errno;
        register int c = getc (stream);

        save_errno = errno;

        /* We always want at least one char left in the buffer, since we
        always (unless we get an error while reading the first char)
        NUL-terminate the line buffer.  */

        assert((*lineptr + *n) == (read_pos + nchars_avail));
        if (nchars_avail < 2) {
            if (*n > min_chunk) {
                *n *= 2;
            } else {
                *n += min_chunk;
            }
            nchars_avail = (int) (*n + *lineptr - read_pos);
            *lineptr = realloc (*lineptr, (size_t)(*n));
            if (!*lineptr) {
                errno = ENOMEM;
                return -1;
            }
            read_pos = *n - nchars_avail + *lineptr;
            assert((*lineptr + *n) == (read_pos + nchars_avail));
        }

        if (ferror (stream) ) {
            /* Might like to return partial line, but there is no
            place for us to store errno.  And we don't want to just
            lose errno.  */
            errno = save_errno;
            return -1;
        }

        if (c == EOF) {
            /* Return partial line, if any.  */
            if (read_pos == *lineptr) {
                return -1;
            } else {
                break;
            }
        }

        *read_pos++ = (char)c;
        nchars_avail--;

        if (c == '\n') {
            /* Return the line.  */
            break;
        }
    }

    /* Done - NUL terminate and return the number of chars read.  */
    *read_pos = '\0';

    ret = (int) (read_pos - *lineptr) ;
    return ret;
}



