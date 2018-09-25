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
int main(int argc, char** argv)
{
    GrB_Index const NUM_NODES = 7;
    GrB_Index const NUM_EDGES = 12;
    GrB_Index row_indices[] = {0, 0, 1, 1, 2, 3, 3, 4, 5, 6, 6, 6};
    GrB_Index col_indices[] = {1, 3, 4, 6, 5, 0, 2, 5, 2, 2, 3, 4};
    bool values[] = {true, true, true, true, true, true,
                     true, true, true, true, true, true};
    GrB_Matrix graph;

    GrB_Matrix_new(&graph, GrB_BOOL, NUM_NODES, NUM_NODES);
    GrB_Matrix_build(graph, row_indices, col_indices, (bool*)values, NUM_EDGES,
                     GrB_LOR);

    pretty_print_matrix_UINT64(graph, "GRAPH");

    // Build a vector to select a source node and another
    // vector to hold the mxv result.
    GrB_Index const SRC_NODE = 0;
    GrB_Vector w, v;
    GrB_Vector_new(&w, GrB_BOOL, NUM_NODES);
    GrB_Vector_new(&v, GrB_BOOL, NUM_NODES);
    GrB_Vector_setElement(w, true, SRC_NODE);

    // Build the transpose/scmp/replace descriptor
    GrB_Descriptor desc;
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_SCMP);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);

    pretty_print_vector_UINT64(w, "wavefront(src)");

    // traverse to neighbors of a frontier iteratively starting with SRC_NODE
    GrB_Index nvals = 0;

    do
    {
        GrB_eWiseAdd(v, GrB_NULL, GrB_NULL,
                     GrB_LOR, v, w, GrB_NULL);
        pretty_print_vector_UINT64(v, "visited");

        GrB_mxv(w, v, GrB_NULL,
                GxB_LOR_LAND_BOOL, graph, w, desc);
        pretty_print_vector_UINT64(w, "wavefront");

        GrB_Vector_nvals(&nvals, w);
    } while (nvals > 0);

    // Cleanup
    GrB_free(&graph);
    GrB_free(&w);
    GrB_free(&v);
    GrB_free(&desc);
    GrB_finalize();
}
