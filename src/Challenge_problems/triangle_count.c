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
 * @file triangle_count.c
 *
 * @brief A triangle counting implementation using GraphBLAS C API.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <GraphBLAS.h>
#include "tutorial_utils.h"

//****************************************************************************
/*
 * Given, L, the lower triangular portion of n x n adjacency matrix A (of and
 * undirected graph), computes the number of triangles in the graph.
 */
uint64_t triangle_count(GrB_Matrix L)             // L: NxN, lower-triangular, bool
{
  GrB_Index n;
  GrB_Matrix_nrows(&n, L);                        // n = # of vertices

  GrB_Matrix C;
  GrB_Matrix_new(&C, GrB_UINT64, n, n);

  GrB_Monoid UInt64Plus;                          // integer plus monoid
  GrB_Monoid_new(&UInt64Plus,GrB_PLUS_UINT64,(uint64_t)0ul);

  GrB_Semiring UInt64Arithmetic;                  // integer arithmetic semiring
  GrB_Semiring_new(&UInt64Arithmetic,UInt64Plus,GrB_TIMES_UINT64);

  GrB_Descriptor desc_tb;                         // Descriptor for mxm
  GrB_Descriptor_new(&desc_tb);
  GrB_Descriptor_set(desc_tb,GrB_INP1,GrB_TRAN); // transpose the second matrix

  GrB_mxm(C, L, GrB_NULL, UInt64Arithmetic, L, L, desc_tb); // C<L> = L *.+ L'

  uint64_t count;
  GrB_reduce(&count, GrB_NULL, UInt64Plus, C, GrB_NULL);    // 1-norm of C

  GrB_free(&C);                      // C matrix no longer needed
  GrB_free(&UInt64Arithmetic);       // Semiring no longer needed
  GrB_free(&UInt64Plus);             // Monoid no longer needed
  GrB_free(&desc_tb);                // descriptor no longer needed

  return count;
}

#if 0
//****************************************************************************
// Logo graph
int main(int argc, char **argv)
{
    GrB_Index const NUM_NODES = 7;
    GrB_Index const NUM_EDGES = 12;
    GrB_Index row_indices[] = {0, 0, 1, 1, 2, 3, 3, 4, 5, 6, 6, 6};
    GrB_Index col_indices[] = {1, 3, 4, 6, 5, 0, 2, 5, 2, 2, 3, 4};

    GrB_Matrix L;
    uint64_t num_triangles = 0;

    // Initialize a GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix_new(&L, GrB_BOOL, NUM_NODES, NUM_NODES);
    // turn the directed graph into an undirected graph and only set lower
    // triangular portion
    for (GrB_Index ix = 0; ix < NUM_EDGES; ++ix)
    {
        if (row_indices[ix] > col_indices[ix])
            GrB_Matrix_setElement(L, true, row_indices[ix], col_indices[ix]);
        else
            GrB_Matrix_setElement(L, true, col_indices[ix], row_indices[ix]);
    }

    num_triangles = triangle_count(L);
    fprintf(stdout, "Number of triangles = %ld\n", (long int)num_triangles);

    return 0;
}
#endif

//****************************************************************************
// Karate club graph
int main(int argc, char **argv)
{
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

    GrB_Matrix L;
    uint64_t num_triangles = 0;

    // Initialize a GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix_new(&L, GrB_BOOL, NUM_NODES, NUM_NODES);
    // turn a directed graph into an undirected graph and only set lower
    // triangular portion
    for (GrB_Index ix = 0; ix < NUM_EDGES; ++ix)
    {
        if (row_indices[ix] > col_indices[ix])
            GrB_Matrix_setElement(L, true, row_indices[ix], col_indices[ix]);
        else
            GrB_Matrix_setElement(L, true, col_indices[ix], row_indices[ix]);
    }

    num_triangles = triangle_count(L);
    fprintf(stdout, "Number of triangles = %ld\n", (long int)num_triangles);

    return 0;
}
