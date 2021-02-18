#ifndef H_XMALLOC_H
#define H_XMALLOC_H
#include <stdlib.h>
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
void *malloc_or_exit(size_t nbytes, const char *file, int line);
#define xmalloc(nbytes) malloc_or_exit((nbytes), __FILE__, __LINE__)
#endif // H_XMALLOC_H
