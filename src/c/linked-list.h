#ifndef H_LINKED_LISTS_H
#define H_LINKED_LISTS_H
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
// implements linked lists.
// 
// a conscell is a data type like 
// a lego block for a snake:
// 
//         .-------.-------.
//         !       !       !
// list -> ! data  ! next  ! -> next conscell
//         !       !       !
//         !-------!-------!
// 
// one part "points" to data,
// other part "points" to next conscell.

/*---------- constructor -----------*/
/*                                  */
typedef struct conscell{
    void *data;
    struct conscell *next;
} conscell;

/*----------- appenders ------------*/
/*                                  */
conscell *ll_push(conscell *list, void *data);
conscell *ll_pop(conscell *list);

/*----------- reorderers -----------*/
/*                                  */
conscell *ll_reverse(conscell *list);

conscell *ll_sort(
    conscell *list,
    int (*cmp) (const void *a, const void *b, void *params),
    void *params);

conscell *ll_filter(conscell *list,
    int (*filter) (const void *a),
    conscell **removed);

/*----------- get length -----------*/
/*                                  */
int ll_length(conscell *list);

/*---------- destructors -----------*/
/*                                  */
void ll_free(conscell *list);

#endif // H_LINKED_LISTS_H
