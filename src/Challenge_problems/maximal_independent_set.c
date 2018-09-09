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
 * @file maximal_independent_set.c
 *
 * @brief A MIS implementation using GraphBLAS C API.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <GraphBLAS.h>
#include "tutorial_utils.h"

// Assign a random number to each element scaled by the inverse of the node's degree.
// This will increase the probability that low degree nodes are selected and larger
// sets are selected.
void setRandom(void *out, const void *in)
{
  uint32_t degree = *(uint32_t*)in;
  double rnd = (double)random()/((double)RAND_MAX);
  *(float*)out = (0.0001f + rnd/(1. + 2.*degree)); // add 1 to prevent divide by zero
}

/*
 * A variant of Luby's randomized algorithm [Luby 1985].
 *
 * Given a numeric n x n adjacency matrix A of an unwieghted and undirected graph (where
 * the value true represents an edge), compute a maximal set of independent vertices and
 * return it in a boolean n-vector, 'iset' where set[i] == true implies vertex i is a member
 * of the set (the iset vector should be uninitialized on input.)
 */
GrB_Info MIS(GrB_Vector *iset, const GrB_Matrix A)
{
  GrB_Index n;
  GrB_Matrix_nrows(&n,A);                       // n = # of rows of A

  GrB_Vector prob;                              // holds random probabilities for each node
  GrB_Vector neighbor_max;                      // holds value of max neighbor probability
  GrB_Vector new_members;                       // holds set of new members to iset
  GrB_Vector new_neighbors;                     // holds set of new neighbors to new iset mbrs.
  GrB_Vector candidates;                        // candidate members to iset

  GrB_Vector_new(&prob,GrB_FP32,n);
  GrB_Vector_new(&neighbor_max,GrB_FP32,n);
  GrB_Vector_new(&new_members,GrB_BOOL,n);
  GrB_Vector_new(&new_neighbors,GrB_BOOL,n);
  GrB_Vector_new(&candidates,GrB_BOOL,n);
  GrB_Vector_new(iset,GrB_BOOL,n);              // Initialize independent set vector, bool

  GrB_Monoid Max;
  GrB_Monoid_new(&Max,GrB_MAX_FP32,0.0f);

  GrB_Semiring maxSelect2nd;                    // Max/Select2nd "semiring"
  GrB_Semiring_new(&maxSelect2nd,Max,GrB_SECOND_FP32);

  GrB_Monoid Lor;
  GrB_Monoid_new(&Lor,GrB_LOR,(bool)false);

  GrB_Semiring Boolean;                         // Boolean semiring
  GrB_Semiring_new(&Boolean,Lor,GrB_LAND);

  // replace
  GrB_Descriptor r_desc;
  GrB_Descriptor_new(&r_desc);
  GrB_Descriptor_set(r_desc,GrB_OUTP,GrB_REPLACE);

  // replace + structural complement of mask
  GrB_Descriptor sr_desc;
  GrB_Descriptor_new(&sr_desc);
  GrB_Descriptor_set(sr_desc,GrB_MASK,GrB_SCMP);
  GrB_Descriptor_set(sr_desc,GrB_OUTP,GrB_REPLACE);

  GrB_UnaryOp set_random;
  GrB_UnaryOp_new(&set_random,setRandom,GrB_FP32,GrB_UINT32);

  // compute the degree of each vertex.
  GrB_Vector degrees;
  GrB_Vector_new(&degrees,GrB_FP64,n);
  GrB_reduce(degrees,GrB_NULL,GrB_NULL,GrB_PLUS_FP64,A,GrB_NULL);

  // Isolated vertices are not candidates: candidates[degrees != 0] = true
  GrB_assign(candidates,degrees,GrB_NULL,true,GrB_ALL,n,GrB_NULL);

  // add all singletons to iset: iset[degree == 0] = 1
  GrB_assign(*iset,degrees,GrB_NULL,true,GrB_ALL,n,sr_desc) ;

  // Iterate while there are candidates to check.
  GrB_Index nvals;
  GrB_Vector_nvals(&nvals, candidates);
  while (nvals > 0) {
    // compute a random probability scaled by inverse of degree
    GrB_apply(prob,candidates,GrB_NULL,set_random,degrees,r_desc);

    // compute the max probability of all neighbors
    GrB_mxv(neighbor_max,candidates,GrB_NULL,maxSelect2nd,A,prob,r_desc);

    // select vertex if its probability is larger than all its active neighbors,
    // and apply a "masked no-op" to remove stored falses
    GrB_eWiseAdd(new_members,GrB_NULL,GrB_NULL,GrB_GT_FP64,prob,neighbor_max,GrB_NULL);
    GrB_apply(new_members,new_members,GrB_NULL,GrB_IDENTITY_BOOL,new_members,r_desc);

    // add new members to independent set.
    GrB_eWiseAdd(*iset,GrB_NULL,GrB_NULL,GrB_LOR,*iset,new_members,GrB_NULL);

    // remove new members from set of candidates c = c & !new
    GrB_eWiseMult(candidates,new_members,GrB_NULL,
                  GrB_LAND,candidates,candidates,sr_desc);

    GrB_Vector_nvals(&nvals, candidates);
    if (nvals == 0) { break; }                  // early exit condition

    // Neighbors of new members can also be removed from candidates
    GrB_mxv(new_neighbors,candidates,GrB_NULL,Boolean,A,new_members,GrB_NULL);
    GrB_eWiseMult(candidates,new_neighbors,GrB_NULL,
                  GrB_LAND,candidates,candidates,sr_desc);

    GrB_Vector_nvals(&nvals, candidates);
  }

  GrB_free(&neighbor_max);                      // free all objects "new'ed"
  GrB_free(&new_members);
  GrB_free(&new_neighbors);
  GrB_free(&prob);
  GrB_free(&candidates);
  GrB_free(&maxSelect2nd);
  GrB_free(&Boolean);
  GrB_free(&Max);
  GrB_free(&Lor);
  GrB_free(&sr_desc);
  GrB_free(&r_desc);
  GrB_free(&set_random);
  GrB_free(&degrees);

  return GrB_SUCCESS;
}

//****************************************************************************
int main(int argc, char **argv)
{
#if 1
    // Logo graph undirected (symmetricized)
    GrB_Index const NUM_NODES = 7;
    GrB_Index const NUM_EDGES = 20;
    GrB_Index row_indices[] = {0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 6};
    GrB_Index col_indices[] = {1, 3, 0, 4, 6, 3, 5, 6, 0, 2, 6, 1, 5, 6, 2, 4, 1, 2, 3, 4};
    bool values[] = {true, true, true, true, true, true, true, true, true, true,
                     true, true, true, true, true, true, true, true, true, true};
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

    bool values[156];
    for (size_t ix = 0; ix < NUM_EDGES; ++ix)
    {
        values[ix] = true;
    }
#endif

    GrB_Matrix graph;

    // Initialize a GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix_new(&graph, GrB_BOOL, NUM_NODES, NUM_NODES);
    GrB_Matrix_build(graph, row_indices, col_indices, (bool*)values, NUM_EDGES,
                     GrB_LOR);

    pretty_print_matrix_FP64(graph, "GRAPH");

    unsigned int const NUM_ITERS = 100;
    GrB_Index set_sizes = 0;
    GrB_Vector iset;
    GrB_Index set_indices[NUM_NODES];
    bool set_vals[NUM_NODES];
    for (unsigned int iter = 0; iter < NUM_ITERS; ++iter)
    {
        srandom(iter); // set a new random seed
        MIS(&iset, graph);
        GrB_Index nvals;
        GrB_Vector_nvals(&nvals, iset);
        set_sizes += nvals;
        fprintf(stdout, "%d: size = %ld:  ", iter, nvals);
        GrB_Vector_extractTuples(set_indices, set_vals, &nvals, iset);
        for (GrB_Index ix = 0; ix < nvals; ++ix)
            fprintf(stdout, " %ld", set_indices[ix]);
        fprintf(stdout, "\n");
        GrB_free(&iset);
    }
    fprintf(stdout, "Average set size = %lf\n", (double)set_sizes/(double)NUM_ITERS);

    // Cleanup
    GrB_free(&graph);
    GrB_finalize();

    return GrB_SUCCESS;
}
