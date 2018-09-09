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
 * @file page_rank.c
 *
 * @brief A PageRank implementation using GraphBLAS C API.  Adapted from the
 *        GraphBLAS Template Library (GBTL) Version 2.0 license.
 */

#include <stdio.h>
#include <GraphBLAS.h>
#include "tutorial_utils.h"

//****************************************************************************
double G_damping_factor;
void times_damping_factor(void *out, void const *in)
{
    *(double*)out = G_damping_factor * (*((double*)in));
}

//****************************************************************************
double G_scaled_teleport;
void plus_scaled_teleport(void *out, void const *in)
{
    *(double*)out = G_scaled_teleport + (*((double*)in));
}

//****************************************************************************
/**
 * @brief Compute the page rank for each node in a graph.
 *
 * Need more documentation
 *
 * @note Because of the random component of this algorithm, the MIS
 *       calculated across various calls to <code>mis</code> may vary.
 *
 * @note This only works with floating point scalars.
 *
 * @param[in]  graph            NxN adjacency matrix (of doubles) of graph to
 *                              compute the page rank.  The structural zero needs
 *                              to be '0' and edges are indicated by '1.0'
 *                              to support use of the Arithmetic FP64 semiring.
 * @param[out] page_rank        N-vector (double) of page ranks (_new already called)
 * @param[in]  damping_factor   The constant to ensure stability in cyclic
 *                              graphs (often 0.85)
 * @param[in]  threshold        The sum of squared errors termination
 *                              threshold.
 * @param[in]  max_iters        The maximum number of iterations to perform (if
 *                              threshold is not met).
 *
 */
GrB_Info page_rank(GrB_Matrix const    graph,
                   GrB_Vector          ranks,
                   double              damping_factor,
                   double              threshold,
                   unsigned int        max_iters)
{
    GrB_Index rows, cols, num_nodes;
    GrB_Matrix_nrows(&rows, graph);
    GrB_Matrix_ncols(&cols, graph);
    GrB_Vector_size(&num_nodes, ranks);

    if ((rows != cols) || (num_nodes != rows))
    {
        return GrB_DIMENSION_MISMATCH;
    }

    // Compute the scaled graph matrix, M, by normalizing the edge
    // weights of the graph by the vertices' out-degrees (num_nodes)
    GrB_Matrix M;
    GrB_Matrix_dup(&M, graph);

    GrB_Vector w;
    GrB_Vector_new(&w, GrB_FP64, num_nodes);
    GrB_reduce(w, GrB_NULL, GrB_NULL, GrB_PLUS_FP64, M, GrB_NULL);
    GrB_apply(w, GrB_NULL, GrB_NULL, GrB_MINV_FP64, w, GrB_NULL);

    GrB_Index num_vals;
    GrB_Vector_nvals(&num_vals, w);

    GrB_Index *indices = malloc(num_vals*sizeof(GrB_Index));
    double *values = malloc(num_vals*sizeof(double));

    GrB_Vector_extractTuples(indices, values, &num_vals, w);

    //populate diagonal:
    GrB_Matrix Mdiag;
    GrB_Matrix_new(&Mdiag, GrB_FP64, num_nodes, num_nodes);
    GrB_Matrix_build(Mdiag, indices, indices, values, num_vals, GrB_PLUS_FP64);

    //Perform matrix multiply to scale rows
    GrB_mxm(M, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_FP64, Mdiag, M, GrB_NULL);

    pretty_print_matrix_FP64(M, "Normalized Graph");

    // scale the normalized edge weights by the damping factor
    G_damping_factor = damping_factor;
    GrB_UnaryOp multiply_damping_factor;
    GrB_UnaryOp_new(&multiply_damping_factor, times_damping_factor,
                    GrB_FP64, GrB_FP64);
    GrB_apply(M, GrB_NULL, GrB_NULL, multiply_damping_factor,
              M, GrB_NULL);

    pretty_print_matrix_FP64(M, "Scaled Graph");
    G_scaled_teleport = (1.0 - damping_factor) / ((double)num_nodes);
    GrB_UnaryOp add_scaled_teleport;
    GrB_UnaryOp_new(&add_scaled_teleport, plus_scaled_teleport,
                    GrB_FP64, GrB_FP64);

    // Initialize all page ranks to 1/N
    GrB_assign(ranks, GrB_NULL, GrB_NULL, 1.0/((double)num_nodes),
               GrB_ALL, num_nodes, GrB_NULL);

    GrB_Vector new_rank, delta;
    GrB_Vector_new(&new_rank, GrB_FP64, num_nodes);
    GrB_Vector_new(&delta, GrB_FP64, num_nodes);
    for (GrB_Index i = 0; i < max_iters; ++i)
    {
        fprintf(stdout, "============= ITERATION %ld ============\n", i);
        pretty_print_vector_FP64(ranks, "rank");

        // Compute the new rank: [1 x N] := [1 x N][N x N]
        GrB_vxm(new_rank, GrB_NULL, GrB_NULL, GxB_PLUS_TIMES_FP64,
                ranks, M, GrB_NULL);
        pretty_print_vector_FP64(new_rank, "step 1:");

        // [1 x N][N x 1] = [1 x 1] = always (1 - damping_factor)
        // rank*(m + scaling_mat*teleport): [1 x 1][1 x N] + [1 x N] = [1 x N]
        // use apply to add "scaled teleport":
        GrB_apply(new_rank, GrB_NULL, GrB_NULL, add_scaled_teleport,
                  new_rank, GrB_NULL);
        pretty_print_vector_FP64(new_rank, "new_rank");

        // Test for convergence - compute squared error
        /// @todo should be mean squared error. (divide r2/N)
        double squared_error;
        GrB_eWiseAdd(delta, GrB_NULL, GrB_NULL, GrB_MINUS_FP64,
                     ranks, new_rank, GrB_NULL);
        GrB_eWiseMult(delta, GrB_NULL, GrB_NULL, GrB_TIMES_FP64,
                      delta, delta, GrB_NULL);
        GrB_reduce(&squared_error, GrB_NULL, GxB_PLUS_FP64_MONOID,
                   delta, GrB_NULL);

        fprintf(stdout, "Squared error = %lf\n", squared_error);

        //copy new page rank vector
        GrB_assign(ranks, GrB_NULL, GrB_NULL, new_rank,
                   GrB_ALL, num_nodes, GrB_NULL);

        // check mean-squared error
        if (squared_error/((double)num_nodes) < threshold)
        {
            break;
        }
    }

    // for any elements missing from page rank vector we need to set
    // to scaled teleport.
    //GrB_assign(new_rank, GrB_NULL, GrB_NULL,
    //           (1.0 - damping_factor) / ((double)num_nodes),
    //           GrB_ALL, GrB_NULL);
    //GrB_eWiseAdd(ranks, ranks, GrB_NULL, GrB_PLUS_FP64,
    //             ranks, new_rank, desc_c);

    free(indices);
    free(values);
    GrB_free(&delta);
    GrB_free(&new_rank);
    GrB_free(&add_scaled_teleport);
    GrB_free(&multiply_damping_factor);
    GrB_free(&M);
    GrB_free(&Mdiag);
    GrB_free(&w);

    return GrB_SUCCESS;
}

//****************************************************************************
int main(int argc, char **argv)
{
#if 1
    // Logo graph
    GrB_Index const NUM_NODES = 7;
    GrB_Index const NUM_EDGES = 12;
    GrB_Index row_indices[] = {0, 0, 1, 1, 2, 3, 3, 4, 5, 6, 6, 6};
    GrB_Index col_indices[] = {1, 3, 4, 6, 5, 0, 2, 5, 2, 2, 3, 4};
    double values[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
#else
    // Karate club graph
    GrB_Index const NUM_NODES = 34;
    GrB_Index const NUM_EDGES = 156;
    GrB_Index row_indices[] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        1,1,1,1,1,1,1,1,1,
        2,2,2,2,2,2,2,2,2,2,
        3,3,3,3,3,3,
        4,4,4,
        5,5,5,5,
        6,6,6,6,
        7,7,7,7,
        8,8,8,8,8,
        9,9,
        10,10,10,
        11,
        12,12,
        13,13,13,13,13,
        14,14,
        15,15,
        16,16,
        17,17,
        18,18,
        19,19,19,
        20,20,
        21,21,
        22,22,
        23,23,23,23,23,
        24,24,24,
        25,25,25,
        26,26,
        27,27,27,27,
        28,28,28,
        29,29,29,29,
        30,30,30,30,
        31,31,31,31,31,31,
        32,32,32,32,32,32,32,32,32,32,32,32,
        33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33};

    GrB_Index col_indices[] = {
        1,2,3,4,5,6,7,8,10,11,12,13,17,19,21,31,
        0,2,3,7,13,17,19,21,30,
        0,1,3,7,8,9,13,27,28,32,
        0,1,2,7,12,13,
        0,6,10,
        0,6,10,16,
        0,4,5,16,
        0,1,2,3,
        0,2,30,32,33,
        2,33,
        0,4,5,
        0,
        0,3,
        0,1,2,3,33,
        32,33,
        32,33,
        5,6,
        0,1,
        32,33,
        0,1,33,
        32,33,
        0,1,
        32,33,
        25,27,29,32,33,
        25,27,31,
        23,24,31,
        29,33,
        2,23,24,33,
        2,31,33,
        23,26,32,33,
        1,8,32,33,
        0,24,25,28,32,33,
        2,8,14,15,18,20,22,23,29,30,31,33,
        8,9,13,14,15,18,19,20,22,23,26,27,28,29,30,31,32};

    double values[156];
    for (size_t ix = 0; ix < NUM_EDGES; ++ix)
    {
        values[ix] = 1.0;
    }
#endif

    GrB_Matrix graph;
    GrB_Vector ranks;

    // Initialize a GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix_new(&graph, GrB_FP64, NUM_NODES, NUM_NODES);
    GrB_Vector_new(&ranks, GrB_FP64, NUM_NODES);
    GrB_Matrix_build(graph, row_indices, col_indices, (double*)values, NUM_EDGES,
                     GrB_LOR);

    pretty_print_matrix_FP64(graph, "GRAPH");
    page_rank(graph, ranks, 0.85, 1.e-5, 100);

    pretty_print_vector_FP64(ranks, "Page Rank");

    // Cleanup
    GrB_free(&ranks);
    GrB_free(&graph);
    GrB_finalize();

    return GrB_SUCCESS;
}
