#include <stdio.h>
#include "xmalloc.h"
// -----------------------------------------------------------------------------
// diego domenzain
// sometime after fall 2017
// 
// code from the book,
// 
// @book{rostamian2014programming,
//   title={Programming Projects in C for Students of Engineering, Science, and Mathematics},
//   author={Rostamian, Rouben},
//   year={2014},
//   publisher={SIAM}
// }
// -----------------------------------------------------------------------------
void *malloc_or_exit(size_t nbytes, const char *file, int line)
{
    void *x;
    if ((x = malloc(nbytes)) == NULL || nbytes == 0)   {
        fprintf(stderr, "%s: line %d: malloc() of %zu bytes failed\n",
                 file, line, nbytes);
        exit(EXIT_FAILURE);
    } else
        return x;
}
