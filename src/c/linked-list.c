#include "linked-list.h"
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

/*----------- appender ------------*/
/*                                 */
conscell *ll_push(conscell *list, void *data)
{
    conscell *new = xmalloc(sizeof *new);
    new -> data = data;     // like new.data = data
    new -> next = list;     // like new.next = list
    return new;
}

/*----------- popper ------------*/
/*                               */
conscell *ll_pop(conscell *list)
{
    if(list != NULL) {
        conscell *p = list->next;
        free(list);
        list = p;
    }
    return list;
}

/*----------- reverse -----------*/
/*                               */
conscell *ll_reverse(conscell *list)
{
    conscell *reverse = NULL;
    while(list){
        reverse = ll_push(reverse, list->data);
        list = ll_pop(list);
    }
    return reverse;
}

/*---------- destructor -----------*/
/*                                 */
void ll_free(conscell *list)
{
    while(list != NULL) {
        conscell *p = list->next;
        free(list);
        list = p;
    }
}
