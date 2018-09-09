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
 * @file matvec.c
 *
 * @brief Code to find sources of a vertex;
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

    // Initialize a GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix graph;
    GrB_Matrix_new(&graph, GrB_BOOL, NUM_NODES, NUM_NODES);
    GrB_Matrix_build(graph, row_indices, col_indices, (bool*)values, NUM_EDGES,
                     GrB_LOR);

    pretty_print_matrix_BOOL(graph, "GRAPH");

    // Build a vector to select a single node in the graph (select)
    // and a vector to hold the result of our operation (result)
    GrB_Index const NODE = 2;
    GrB_Vector select, result;
    GrB_Vector_new(&select, GrB_BOOL, NUM_NODES);
    GrB_Vector_new(&result, GrB_BOOL, NUM_NODES);
    GrB_Vector_setElement(select, true, NODE);

    // find source vertices to NODE

    pretty_print_vector_BOOL(select, "Target node");
    GrB_mxv(result, GrB_NULL, GrB_NULL,
            GxB_LOR_LAND_BOOL, graph, select, GrB_NULL);
    pretty_print_vector_BOOL(result, "sources");

    // Check results
    {
        bool error_found = false;
        GrB_Index nvals;
        GrB_Vector_nvals(&nvals, result);
        if (nvals != 3)
        {
            fprintf(stderr, "ERROR: wrong number of sources (!= 2): %ld\n",
                    (long)nvals);
            error_found = true;
        }

        bool val;

        if (GrB_Vector_extractElement(&val, result, 3UL))
        {
            fprintf(stderr, "ERROR: missing source 3.\n");
            error_found = true;
        }
        if (GrB_Vector_extractElement(&val, result, 5UL))
        {
            fprintf(stderr, "ERROR: missing source 5.\n");
            error_found = true;
        }
        if (GrB_Vector_extractElement(&val, result, 6UL))
        {
            fprintf(stderr, "ERROR: missing source 6.\n");
            error_found = true;
        }

        if (!error_found)
        {
            fprintf(stderr, "GrB_mxv test passed.\n");
        }
    }

    // Cleanup
    GrB_free(&graph);
    GrB_free(&select);
    GrB_free(&result);
    GrB_finalize();
}
