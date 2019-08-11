/*
    Copyright 2019 xxx.

    All Rights Reserved.

    NO WARRANTY. THIS MATERIAL IS FURNISHED ON AN "AS-IS" BASIS. THE LAGRAPH
    CONTRIBUTORS MAKE NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
    AS TO ANY MATTER INCLUDING, BUT NOT LIMITED TO, WARRANTY OF FITNESS FOR
    PURPOSE OR MERCHANTABILITY, EXCLUSIVITY, OR RESULTS OBTAINED FROM USE OF
    THE MATERIAL. THE CONTRIBUTORS DO NOT MAKE ANY WARRANTY OF ANY KIND WITH
    RESPECT TO FREEDOM FROM PATENT, TRADEMARK, OR COPYRIGHT INFRINGEMENT.

    Released under a BSD license, please see the LICENSE file distributed with
    this Software or contact permission@sei.cmu.edu for full terms.

    Created, in part, with funding and support from the United States
    Government.  (see Acknowledgments.txt file).

    This program includes and/or can make use of certain third party source
    code, object code, documentation and other files ("Third Party Software").
    See LICENSE file for more details.

*/

//------------------------------------------------------------------------------

// Contributed by Scott McMillan/CMU, Tim Mattson/Intel

#include <stdio.h>
#include "LAGraph.h"
#include "tutorial_utils.h"

#define OK(method)                                                          \
{                                                                           \
    GrB_Info this_info = method ;                                           \
    if (! (this_info == GrB_SUCCESS || this_info == GrB_NO_VALUE))          \
    {                                                                       \
        printf ("GrB failure: [%d] %s\n", this_info, GrB_error ( )) ;       \
        return (this_info) ;                                                \
    }                                                                       \
}

//------------------------------------------------------------------------------
// read_mtx main:
//------------------------------------------------------------------------------
int main (int argc, char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Error\nUsage: %s <matrix_market_file.mtx>\n", argv[0]);
        return 1;
    }

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS )
    printf ("testing LAGraph_xinit (requires SuiteSparse:GraphBLAS)\n") ;
    LAGraph_xinit (malloc, calloc, realloc, free, true) ;
    #else
    printf ("LAGraph_init\n") ;
    LAGraph_init ( ) ;
    #endif

    printf("reading input graph: %s\n", argv[1]);
    GrB_Matrix A = NULL;
    FILE *fd = fopen(argv[1], "r");

    double tic[2], t;
    LAGraph_tic (tic);
    OK (LAGraph_mmread (&A, fd); )
    t = LAGraph_toc(tic) ;
    if (fd != NULL) fclose(fd);
    printf ("time taken: %g sec\n", t) ;

    GrB_Index num_rows, num_cols, num_vals;
    GrB_Matrix_nrows(&num_rows, A);
    GrB_Matrix_ncols(&num_cols, A);
    GrB_Matrix_nvals(&num_vals, A);
    printf ("nrows/ncols/nvals = %ld/%ld/%ld\n", num_rows, num_cols, num_vals);

    GrB_Vector degree;
    GrB_Vector_new(&degree, GrB_UINT64, num_rows);
    GrB_reduce(degree, GrB_NULL, GrB_NULL, GrB_PLUS_UINT64, A, GrB_NULL);

    GrB_Index max_val, max_index;
    for (GrB_Index ix = 0; ix < num_rows; ++ix)
    {
        GrB_Index val;
        if (GrB_NO_VALUE != GrB_Vector_extractElement(&val, degree, ix))
        {
            if (val > max_val)
            {
                max_val   = val;
                max_index = ix;
            }
        }
    }

    printf("Author with the most coauthor/paper combos: %ld (count: %ld)\n", max_index, max_val);

    // Use LACC connected components algorithm
    GrB_Vector LACC_result;
    LAGraph_lacc(A, &LACC_result);
    pretty_print_vector_UINT64(LACC_result, "Connected components");

    GrB_Index max_component = 666;
    GrB_Vector_extractElement(&max_component, LACC_result, max_index);
    printf("Component for author %ld: %ld\n", max_index, max_component);

    GrB_free(&A);
    LAGraph_finalize ( ) ;
    return (GrB_SUCCESS) ;
}
