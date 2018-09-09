/*
 * This file is part of the GraphBLAS Tutorial materials,
 * Copyright (c) 2018 Carnegie Mellon University and Intel Corporation.
 * All Rights Reserved
 *
 * THIS SOFTWARE IS PROVIDED "AS IS," WITH NO WARRANTIES WHATSOEVER. CARNEGIE
 * MELLON UNIVERSITY AND INTEL CORPORATION EXPRESSLY DISCLAIMS TO THE FULLEST
 * EXTENT PERMITTED BY LAW ALL EXPRESS, IMPLIED, AND STATUTORY WARRANTIES,
 * INCLUDING, WITHOUT LIMITATION, THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF PROPRIETARY RIGHTS.
 *
 * Released under a BSD (SEI)-style license, please see LICENSE.txt for
 * full terms.
 *
 * DM18-xxx
 *
 * Authors: Scott McMillan, Timothy G. Mattson
 */

/**
 * @file bfs_level.c
 *
 * @brief A level BFS implementation using GraphBLAS C API.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <GraphBLAS.h>
#include "tutorial_utils.h"

/*
 * Given a boolean n x n adjacency matrix, A, and a source vertex, src,
 * performs a BFS traversal of the graph and sets levels[i] to the level in
 * which vertex i is visited (levels[src] == 1).  If i is not reacheable from
 * s, then levels[i] = 0. (Vector v should be empty/uninitialized on input.)
 */
GrB_Info bfs_level(GrB_Vector *levels, GrB_Matrix A, GrB_Index src)
{
    GrB_Index n;
    GrB_Matrix_nrows(&n, A);                      // n = # of rows of A

    GrB_Vector_new(levels, GrB_INT32, n);         // Vector<int32_t> levels(n)

    GrB_Vector w;                                 // wavefront vertices visited in each level
    GrB_Vector_new(&w, GrB_BOOL, n);              // Vector<bool> w(n)
    GrB_Vector_setElement(w, (bool)true, src);    // w[src] = true, false elsewhere

    GrB_Monoid Lor;                               // Logical-or monoid
    GrB_Monoid_new(&Lor, GrB_LOR, (bool)false);

    GrB_Semiring Boolean;                         // Boolean semiring
    GrB_Semiring_new(&Boolean, Lor, GrB_LAND);

    GrB_Descriptor desc;                          // Descriptor for vxm: replace+scmp
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_SCMP);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);

    /*
     * BFS traversal and label the vertices.
     */
    int32_t level = 0;                            // level in BFS traversal
    GrB_Index nvals = 0;                          // nvals == 0 when no successor found
    do {
        ++level;                                  // next level (start with 1)
        GrB_assign(*levels, w, GrB_NULL,
                   level, GrB_ALL, n, GrB_NULL);  // levels[w] = d
        GrB_vxm(w, *levels, GrB_NULL, Boolean,
                w, A, desc);                      // w[!levels] = w ||.&& A ; finds all the

        GrB_Vector_nvals(&nvals, w);
    } while (nvals);                              // if there is no successor in w, we are done.

    GrB_free(&w);                                 // Cleanup
    GrB_free(&Lor);
    GrB_free(&Boolean);
    GrB_free(&desc);

    return GrB_SUCCESS;
}

//****************************************************************************
// Logo graph
int main(int argc, char **argv)
{
    GrB_Index const NUM_NODES = 7;
    GrB_Index const NUM_EDGES = 12;
    GrB_Index row_indices[] = {0, 0, 1, 1, 2, 3, 3, 4, 5, 6, 6, 6};
    GrB_Index col_indices[] = {1, 3, 4, 6, 5, 0, 2, 5, 2, 2, 3, 4};
    bool values[] = {true, true, true, true, true, true,
                     true, true, true, true, true, true};
    GrB_Matrix graph;
    GrB_Vector levels;

    // Initialize a GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix_new(&graph, GrB_BOOL, NUM_NODES, NUM_NODES);
    GrB_Matrix_build(graph, row_indices, col_indices, (bool*)values, NUM_EDGES,
                     GrB_LOR);

    pretty_print_matrix_BOOL(graph, "GRAPH");

    GrB_Index const SRC_NODE = 1;
    bfs_level(&levels, graph, SRC_NODE);

    pretty_print_vector_FP64(levels, "level vector (src == 1)");

    // Cleanup
    GrB_free(&levels);
    GrB_free(&graph);
    GrB_finalize();

    return GrB_SUCCESS;
}
