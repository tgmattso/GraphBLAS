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

// Run the following from the src directory:
//
// ./AnalyzeGraph.exe Data/hpec_coauthors.mtx

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

// From LACC_GraphBLAS.c in LAGraph
GrB_Info CountCC(GrB_Vector parents, GrB_Index* countcc);

//------------------------------------------------------------------------------
GrB_Index G_target_id = 0;

void eq_target(void *z, const void *x)
{
    (*(bool*)z) = (bool)((*(GrB_Index*)x) == G_target_id);
    //printf("target %ld : element %ld => %d\n", G_target_id, (*(GrB_Index*)x),
    //       (*(bool*)z));
}

//------------------------------------------------------------------------------
// BFS visited: (from matvecTransIterVisitedExitFlag.c)
//------------------------------------------------------------------------------
GrB_Info BFS_mark(GrB_Matrix const A, GrB_Index src_node, GrB_Vector v)
{
    GrB_Index n;
    GrB_Matrix_nrows(&n, A);
    GrB_Vector w;
    GrB_Vector_new(&w, GrB_BOOL, n);  // wavefront
    GrB_Vector_setElement(w, true, src_node);

    GrB_Descriptor desc;                          // Descriptor for vxm: replace+scmp
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_SCMP);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);

    // traverse to neighbors of a frontier iteratively starting with SRC_NODE
    GrB_Index nvals = 0;

    do
    {
        GrB_eWiseAdd(v, GrB_NULL, GrB_NULL, GrB_LOR, v, w, GrB_NULL);
        GrB_vxm(w, v, GrB_NULL, GxB_LOR_LAND_BOOL, w, A, desc);
        GrB_Vector_nvals(&nvals, w);
    } while (nvals > 0);

    return GrB_SUCCESS;
}

//------------------------------------------------------------------------------
// connected_components:
//------------------------------------------------------------------------------
GrB_Info CC_iterative(GrB_Matrix const  A,
                      GrB_Vector       *components,
                      GrB_Index        *num_components)
{
    GrB_Index n;
    GrB_Matrix_nrows(&n, A);
    GrB_Vector_new(components, GrB_UINT64, n);

    GrB_Vector v;
    GrB_Vector_new(&v, GrB_BOOL, n);  // visited list

    *num_components = 0;
    GrB_Index max_component_num = 0;
    GrB_Index max_component_size = 0;

    for (GrB_Index src_node = 0; src_node < n; ++src_node)
    {
        GrB_Index c_num;
        if (GrB_NO_VALUE == GrB_Vector_extractElement(&c_num, *components, src_node))
        {
            // Traverse as far as you can go from src_node marking all visited nodes
            BFS_mark(A, src_node, v);

            // Merge visited list into components (give component source node name)
            // components[v] = src_node
            GrB_assign(*components, v, GrB_NULL, src_node, GrB_ALL, n, GrB_NULL);
            ++(*num_components);

            // Just for curiousity...not needed...
            GrB_Index component_size;
            GrB_Vector_nvals(&component_size, v);
            //printf("Component %ld: num nodes = %ld\n", src_node, component_size);
            if (component_size > max_component_size)
            {
                max_component_size = component_size;
                max_component_num = src_node;
            }

            GrB_Vector_clear(v);
        }
    }

    printf("Largest component #%ld (size = %ld)\n", max_component_num, max_component_size);
    GrB_free(&v);
    return GrB_SUCCESS;
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

    GrB_Index max_val, target_index;
    for (GrB_Index ix = 0; ix < num_rows; ++ix)
    {
        GrB_Index val;
        if (GrB_NO_VALUE != GrB_Vector_extractElement(&val, degree, ix))
        {
            if (val > max_val)
            {
                max_val   = val;
                target_index = ix;
            }
        }
    }

    printf("Author with the most coauthor/paper combos: %ld (count: %ld)\n", target_index, max_val);

    //============================================================
    // Use LAGraph's LACC connected components algorithm
    printf("#1: Running LAGraphs LACC algorithm\n");
    GrB_Vector LACC_result;
    LAGraph_tic (tic);
    LAGraph_lacc(A, &LACC_result);
    t = LAGraph_toc(tic) ;
    printf ("time taken: %g sec\n", t) ;

    // get the number of components found in LACC
    GrB_Index num_components;
    CountCC(LACC_result, &num_components);
    printf("Number of connected components: %ld\n", num_components);

    //pretty_print_vector_UINT64(LACC_result, "Connected components #1");

    GrB_Index cluster_num = 666;
    GrB_Vector_extractElement(&cluster_num, LACC_result, target_index);
    printf("Component for author %ld: %ld\n", target_index, cluster_num);

    //============================================================
    // Run the class generated CC algorithm above
    printf("#2: Running naive CC algorithm\n");
    GrB_Vector components;
    GrB_Index  num_comps;
    LAGraph_tic (tic);
    CC_iterative(A, &components, &num_comps);
    t = LAGraph_toc(tic) ;
    printf ("time taken: %g sec\n", t) ;

    printf("Number of connected components: %ld\n", num_comps);

    //pretty_print_vector_UINT64(components, "Connected components #2");

    GrB_Vector_extractElement(&cluster_num, components, target_index);
    printf("Component for author %ld: %ld\n", target_index, cluster_num);

    //===========================================================
    // Find all the elements in the target cluster
    G_target_id = cluster_num;
    GrB_UnaryOp target_eq = NULL;
    GrB_UnaryOp_new(&target_eq, &eq_target, GrB_BOOL, GrB_UINT64);

    GrB_Vector cluster_mask;
    GrB_Vector_new(&cluster_mask, GrB_BOOL, num_rows);

    // Find all elements belonging to cluster labeled cluster_num
    GrB_apply(cluster_mask, GrB_NULL, GrB_NULL,
              target_eq, components, GrB_NULL);

    GrB_Index nc;
    GrB_Vector_nvals(&nc, cluster_mask);
    printf("Cluster mask nvals (before masking): %ld\n", nc);

    // remove all elements that have 'stored false'
    GrB_Descriptor desc;
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);

    GrB_apply(cluster_mask, cluster_mask, GrB_NULL,
              GrB_IDENTITY_BOOL, cluster_mask, desc);

    GrB_Vector_nvals(&nc, cluster_mask);
    printf("Cluster mask nvals (after masking): %ld\n", nc);

    // get number of elements in the mask and extract the indices
    GrB_Index component_size;
    GrB_Vector_nvals(&component_size, cluster_mask);
    printf("Component size: %ld\n", component_size);
    GrB_Index *cluster_indices = malloc(component_size*sizeof(GrB_Index));
    bool      *vals = malloc(component_size*sizeof(bool));

    GrB_Vector_extractTuples(cluster_indices, vals,
                             &component_size, cluster_mask);
    printf("Extracted component size: %ld\n", component_size);

    //===========================================================
    // extract the component from the rest of the graph
    GrB_Matrix A_comp;
    GrB_Matrix_new(&A_comp, GrB_UINT64, component_size, component_size);
    GrB_extract(A_comp, GrB_NULL, GrB_NULL, A,
                cluster_indices, component_size,
                cluster_indices, component_size,
                GrB_NULL);

    GrB_Vector pr;
    double max_rank = 0.0;
    GrB_Index max_rank_id;
    LAGraph_pagerank2(&pr, A_comp, 0.85, 100);
    printf("ID:\tRank\n");
    for (GrB_Index ix = 0; ix < component_size; ++ix)
    {
        double rank=-1.;
        GrB_Vector_extractElement(&rank, pr, ix);
        printf("%ld\t%lf\n", cluster_indices[ix], rank);
        if (rank > max_rank)
        {
            max_rank = rank;
            max_rank_id = cluster_indices[ix];
        }
    }

    printf("Author with the highest rank: %ld (%lf)\n", max_rank_id, max_rank);

    //===========================================================
    free(cluster_indices);
    free(vals);
    GrB_free(&cluster_mask);
    GrB_free(&target_eq);
    GrB_free(&LACC_result);
    GrB_free(&components);
    GrB_free(&A);
    LAGraph_finalize ( ) ;
    return (GrB_SUCCESS) ;
}
