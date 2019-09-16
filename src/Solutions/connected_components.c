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
 * @file matvecTransIterExitFlag_v2.c
 *
 * @brief Code to hop to neighbors of a set of vertices iteratively and exit
 *        when the frontier becomes empty.
 *
 */

#include <stdio.h>
#include <GraphBLAS.h>
#include "tutorial_utils.h"

//****************************************************************************
GrB_Info BFS(GrB_Matrix const graph,
             GrB_Index        src_node,
             GrB_Vector       v)
{
    GrB_Index num_nodes;
    GrB_Matrix_nrows(&num_nodes, graph);
    GrB_Vector w;
    GrB_Vector_new(&w, GrB_BOOL, num_nodes);
    GrB_Vector_setElement(w, true, src_node);

    GrB_Descriptor desc;    // Descriptor for vxm: replace+scmp+trans
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_SCMP);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);

    pretty_print_vector_UINT64(w, "wavefront(src)");

    // traverse to neighbors of a frontier iteratively starting with src_node
    GrB_Index nvals = 0;

    do
    {
        GrB_eWiseAdd(v, GrB_NULL, GrB_NULL, GrB_LOR, v, w, GrB_NULL);
        pretty_print_vector_UINT64(v, "visited");

        GrB_mxv(w, v, GrB_NULL, GxB_LOR_LAND_BOOL, graph, w, desc);
        pretty_print_vector_UINT64(w, "wavefront");

        GrB_Vector_nvals(&nvals, w);
    } while (nvals > 0);

    GrB_free(&w);
    GrB_free(&desc);
    return GrB_SUCCESS;
}

//****************************************************************************
int main(int argc, char** argv)
{
    GrB_init(GrB_BLOCKING);

    GrB_Index const NUM_NODES = 14;
    GrB_Index const NUM_EDGES = 40;
    GrB_Index row_indices[] =
        {0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 6,
         7, 7, 8, 8, 8, 9, 9, 9,10,10,10,11,11,11,12,12,13,13,13,13};
    GrB_Index col_indices[] =
        {1, 3, 0, 4, 6, 3, 5, 6, 0, 2, 6, 1, 5, 6, 2, 4, 1, 2, 3, 4,
         8,10, 7,11,13,10,12,13, 7, 9,13, 8,12,13, 9,11, 8, 9,10,11};

    bool values[] = {true, true, true, true, true, true, true, true,
                     true, true, true, true, true, true, true, true,
                     true, true, true, true, true, true, true, true,
                     true, true, true, true, true, true, true, true,
                     true, true, true, true, true, true, true, true};
    GrB_Matrix graph;

    GrB_Matrix_new(&graph, GrB_BOOL, NUM_NODES, NUM_NODES);
    GrB_Matrix_build(graph, row_indices, col_indices, (bool*)values,
                     NUM_EDGES, GrB_LOR);

    pretty_print_matrix_UINT64(graph, "GRAPH");

    // ------------------ connected components algorithm ------------------
    GrB_Index tmp = 0, num_ccs = 0;
    GrB_Vector cc_ids, visited;
    GrB_Vector_new(&cc_ids, GrB_UINT64, NUM_NODES);
    GrB_Vector_new(&visited, GrB_BOOL, NUM_NODES);

    // Build a vector to select a source node and another
    // vector to hold the mxv result.
    for (GrB_Index src = 0; src < NUM_NODES; ++src)
    {
        if (GrB_NO_VALUE ==
            GrB_Vector_extractElement(&tmp, cc_ids, src))
        {
            // Traverse from src_node marking all visited nodes
            BFS(graph, src, visited);

            // Merge visited list into cc_ids (use source node as ID)
            // cc_ids[visited] = src
            GrB_assign(cc_ids, visited, GrB_NULL,
                       src, GrB_ALL, NUM_NODES, GrB_NULL);
            ++num_ccs;
        }
        GrB_Vector_clear(visited);
    }

    printf("Number of connected components: %ld\n", (unsigned long)num_ccs);
    pretty_print_vector_UINT64(cc_ids, "CC IDs");

    // Cleanup
    GrB_free(&graph);
    GrB_free(&visited);
    GrB_free(&cc_ids);
    GrB_finalize();
}
