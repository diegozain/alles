#ifndef H_ARRAY_H
#define H_ARRAY_H
#include <stdlib.h>
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
/*
 *                        _
 *      vector           | |
 *                       | | n
 *                       | |
 *                       |_| 
 */
// -----------------------------------------------------------------------------
// vector maker and destructor
#define make_vector(v,n) ((v) = xmalloc((n) * sizeof *(v)))
#define free_vector(v) do { free(v); v = NULL; } while (0)
// printer vector
#define print_vector(fmt, v, n) do {\
    size_t print_vector_loop_counter;\
    for (print_vector_loop_counter = 0;\
            print_vector_loop_counter < (n);\
            print_vector_loop_counter++)\
        printf(fmt, (v)[print_vector_loop_counter]);\
    putchar('\n');\
} while (0)

// -----------------------------------------------------------------------------
/*                         
 *                        _n________
 *      matrix           |          |
 *                       |          | m
 *                       |          |
 *                       |__________| 
 */
// -----------------------------------------------------------------------------
// matrix maker (m x n) and destructor
#define make_matrix(a, m, n) do {\
    size_t make_matrix_loop_counter;\
    make_vector(a, (m) + 1);\
    for (make_matrix_loop_counter = 0;\
            make_matrix_loop_counter < (m);\
            make_matrix_loop_counter++)\
        make_vector((a)[make_matrix_loop_counter], (n));\
    (a)[m] = NULL;\
} while (0)
// free the matrix
#define free_matrix(a) do {\
    if (a != NULL) {\
        size_t make_matrix_loop_counter;\
        for (make_matrix_loop_counter = 0;\
                (a)[make_matrix_loop_counter] != NULL;\
                make_matrix_loop_counter++)\
            free_vector((a)[make_matrix_loop_counter]);\
        free_vector(a);\
    }\
} while (0)
// printer matrix
#define print_matrix(fmt, a, m, n) do {\
    size_t print_matrix_loop_counter;\
    for (print_matrix_loop_counter = 0;\
            print_matrix_loop_counter < (m);\
            print_matrix_loop_counter++)\
        print_vector(fmt, a[print_matrix_loop_counter], n);\
    putchar('\n');\
} while (0)
// -----------------------------------------------------------------------------
/*                         ___________
 *                        /          /|
 *                       /__________/ |
 *      cube             |          | | l
 *                       |          | |
 *                       |          | /
 *                       |__________|/ m
 *                             n           
 */                         
// -----------------------------------------------------------------------------
// cube maker (m x n x l) and destructor
#define make_cube(a, m, n, l) do {\
    size_t make_cube_loop_counter;\
    make_vector(a, (m) + 1);\
    for (make_cube_loop_counter = 0;\
            make_cube_loop_counter < (m);\
            make_cube_loop_counter++)\
        make_matrix((a)[make_cube_loop_counter], (n), (l));\
    (a)[m] = NULL;\
} while (0)
// free the cube
#define free_cube(a) do {\
    if (a != NULL) {\
        size_t make_cube_loop_counter;\
        for (make_cube_loop_counter = 0;\
                (a)[make_cube_loop_counter] != NULL;\
                make_cube_loop_counter++)\
            free_matrix((a)[make_cube_loop_counter]);\
        free_vector(a);\
    }\
} while (0)
// printer cube
#define print_cube(fmt, a, m, n, l) do {\
    size_t print_cube_counter;\
    for (print_cube_counter = 0;\
            print_cube_counter < (m);\
            print_cube_counter++)\
        print_matrix(fmt, a[print_cube_counter], n, l);\
    putchar('\n');\
} while (0)
// -----------------------------------------------------------------------------
/*                         ___________            ___________
 *                        /          /|          /          /| 
 *                       /__________/ |         /__________/ |
 *  hyper-cube           |          | | p       |          | | p  
 *                       |          | |   ...   |          | | 
 *                       |          | /         |          | / 
 *                       |__________|/ n        |__________|/ n
 *                             l                      l
 *                           ------------------------
 *                                      m
 */                         
// -----------------------------------------------------------------------------
// hypcube maker (m x n x l) and destructor
#define make_hypcube(a, m, n, l, p) do {\
    size_t make_hypcube_loop_counter;\
    make_vector(a, (m) + 1);\
    for (make_hypcube_loop_counter = 0;\
            make_hypcube_loop_counter < (m);\
            make_hypcube_loop_counter++)\
        make_cube((a)[make_hypcube_loop_counter], (n), (l), (p));\
    (a)[m] = NULL;\
} while (0)
// free the hypcube
#define free_hypcube(a) do {\
    if (a != NULL) {\
        size_t make_hypcube_loop_counter;\
        for (make_hypcube_loop_counter = 0;\
                (a)[make_hypcube_loop_counter] != NULL;\
                make_hypcube_loop_counter++)\
            free_cube((a)[make_hypcube_loop_counter]);\
        free_vector(a);\
    }\
} while (0)
// printer hypcube
#define print_hypcube(fmt, a, m, n, l, p) do {\
    size_t print_hypcube_counter;\
    for (print_hypcube_counter = 0;\
            print_hypcube_counter < (m);\
            print_hypcube_counter++)\
        print_cube(fmt, a[print_hypcube_counter], n, l, p);\
    putchar('\n');\
} while (0)
// -----------------------------------------------------------------------------
#endif // H_ARRAY_H
