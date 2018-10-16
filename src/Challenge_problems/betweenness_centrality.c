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
 * @file betweenness_centrality.c
 *
 * @brief Two implementations of the betweenness centrality algorithm using
 *        the GraphBLAS C API. One for single source contribution and one for
 *        batch contribution.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <GraphBLAS.h>
#include "tutorial_utils.h"

//****************************************************************************
/*
 * Given a boolean n x n adjacency matrix A and a source vertex s,
 * compute the BC-metric vector delta, which should be empty on input.
 */
GrB_Info BC(GrB_Vector *delta, GrB_Matrix A, GrB_Index s)
{
  GrB_Index n;
  GrB_Matrix_nrows(&n,A);                       // n = # of vertices in graph

  GrB_Vector_new(delta,GrB_FP32,n);             // Vector<float> delta(n)

  GrB_Matrix sigma;                             // Matrix<int32_t> sigma(n,n)
  GrB_Matrix_new(&sigma,GrB_INT32,n,n);         // sigma[d,k] = #shortest paths to node k at level d

  GrB_Vector q;
  GrB_Vector_new(&q, GrB_INT32, n);             // Vector<int32_t> q(n) of path counts
  GrB_Vector_setElement(q,1,s);                 // q[s] = 1

  GrB_Vector p;                                 // Vector<int32_t> p(n) shortest path counts so far
  GrB_Vector_dup(&p, q);                        // p = q

  GrB_Monoid Int32Add;                          // Monoid <int32_t,+,0>
  GrB_Monoid_new(&Int32Add,GrB_PLUS_INT32,0);
  GrB_Semiring Int32AddMul;                     // Semiring <int32_t,int32_t,int32_t,+,*,0,1>
  GrB_Semiring_new(&Int32AddMul,Int32Add,GrB_TIMES_INT32);

  GrB_Descriptor desc;                          // Descriptor for vxm
  GrB_Descriptor_new(&desc);
  GrB_Descriptor_set(desc,GrB_MASK,GrB_SCMP);   // structural complement of the mask
  GrB_Descriptor_set(desc,GrB_OUTP,GrB_REPLACE);// clear the output before assignment

  GrB_Descriptor tr1;                           // Transpose 1st input argument
  GrB_Descriptor_new(&tr1);
  GrB_Descriptor_set(tr1,GrB_INP0,GrB_TRAN);    // structural complement of the mask

  GrB_vxm(q,p,GrB_NULL,Int32AddMul,q,A,desc);   // get the first set of out neighbors

  /*
   * BFS phase
   */
  int32_t d = 0;                                // BFS level number
  int32_t sum = 0;                              // sum == 0 when BFS phase is complete
  do {
    GrB_assign(sigma,GrB_NULL,GrB_NULL,q,d,GrB_ALL,n,GrB_NULL); // sigma[d,:] = q
    GrB_eWiseAdd(p,GrB_NULL,GrB_NULL,Int32AddMul,p,q,GrB_NULL); // accumulate path counts on this level
    GrB_vxm(q,p,GrB_NULL,Int32AddMul,q,A,desc);                 // q = # paths to nodes reachable
                                                                //    from current level
    GrB_reduce(&sum,GrB_NULL,Int32Add,q,GrB_NULL);              // sum path counts at this level
    ++d;
  } while (sum);

  /*
   * BC computation phase
   * (t1,t2,t3,t4) are temporary vectors
   */
  GrB_Monoid FP32Add;                           // Monoid <float,float,float,+,0.0>
  GrB_Monoid_new(&FP32Add,GrB_PLUS_FP32,0.0f);

  GrB_Monoid FP32Mul;                           // Monoid <float,float,float,*,1.0>
  GrB_Monoid_new(&FP32Mul,GrB_TIMES_FP32,1.0f);

  GrB_Semiring FP32AddMul;                      // Semiring <float,float,float,+,*,0.0,1.0>
  GrB_Semiring_new(&FP32AddMul,FP32Add,GrB_TIMES_FP32);

  GrB_Vector t1; GrB_Vector_new(&t1,GrB_FP32,n);
  GrB_Vector t2; GrB_Vector_new(&t2,GrB_FP32,n);
  GrB_Vector t3; GrB_Vector_new(&t3,GrB_FP32,n);
  GrB_Vector t4; GrB_Vector_new(&t4,GrB_FP32,n);
  for(int i=d-1; i>0; i--)
  {
    GrB_assign(t1,GrB_NULL,GrB_NULL,1.0f,GrB_ALL,n,GrB_NULL);          // t1 = 1+delta
    GrB_eWiseAdd(t1,GrB_NULL,GrB_NULL,FP32Add,t1,*delta,GrB_NULL);
    GrB_extract(t2,GrB_NULL,GrB_NULL,sigma,GrB_ALL,n,i,tr1);           // t2 = sigma[i,:]
    GrB_eWiseMult(t2,GrB_NULL,GrB_NULL,GrB_DIV_FP32,t1,t2,GrB_NULL);   // t2 = (1+delta)/sigma[i,:]
    GrB_mxv(t3,GrB_NULL,GrB_NULL,FP32AddMul,A,t2,GrB_NULL);            // add contributions made by
                                                                       //     successors of a node
    GrB_extract(t4,GrB_NULL,GrB_NULL,sigma,GrB_ALL,n,i-1,tr1);         // t4 = sigma[i-1,:]
    GrB_eWiseMult(t4,GrB_NULL,GrB_NULL,FP32Mul,t4,t3,GrB_NULL);        // t4 = sigma[i-1,:]*t3
    GrB_eWiseAdd(*delta,GrB_NULL,GrB_NULL,FP32Add,*delta,t4,GrB_NULL); // accumulate into delta
  }

  GrB_free(&sigma);
  GrB_free(&q); GrB_free(&p);
  GrB_free(&Int32AddMul); GrB_free(&Int32Add); GrB_free(&FP32AddMul);
  GrB_free(&FP32Add); GrB_free(&FP32Mul);
  GrB_free(&desc);
  GrB_free(&t1); GrB_free(&t2); GrB_free(&t3); GrB_free(&t4);

  return GrB_SUCCESS;
}

//****************************************************************************
// Compute partial BC metric for a subset of source vertices, s, in graph A
GrB_Info BC_update(GrB_Vector *delta, GrB_Matrix A, GrB_Index *s, GrB_Index nsver)
{
  GrB_Index n;
  GrB_Matrix_nrows(&n, A);                             // n = # of vertices in graph
  GrB_Vector_new(delta,GrB_FP32,n);                    // Vector<float> delta(n)

  GrB_Monoid Int32Add;                                 // Monoid <int32_t,+,0>
  GrB_Monoid_new(&Int32Add,GrB_PLUS_INT32,0);
  GrB_Semiring Int32AddMul;                            // Semiring <int32_t,int32_t,int32_t,+,*,0>
  GrB_Semiring_new(&Int32AddMul,Int32Add,GrB_TIMES_INT32);

  // Descriptor for BFS phase mxm
  GrB_Descriptor desc_tsr;
  GrB_Descriptor_new(&desc_tsr);
  GrB_Descriptor_set(desc_tsr,GrB_INP0,GrB_TRAN);      // transpose the adjacency matrix
  GrB_Descriptor_set(desc_tsr,GrB_MASK,GrB_SCMP);      // complement the mask
  GrB_Descriptor_set(desc_tsr,GrB_OUTP,GrB_REPLACE);   // clear output before result is stored

  // index and value arrays needed to build numsp
  GrB_Index *i_nsver = (GrB_Index*)malloc(sizeof(GrB_Index)*nsver);
  int32_t   *ones    = (int32_t*)  malloc(sizeof(int32_t)*nsver);
  for(int i=0; i<nsver; ++i) {
    i_nsver[i] = i;
    ones[i] = 1;
  }

  // numsp: structure holds the number of shortest paths for each node and starting vertex
  // discovered so far.  Initialized to source vertices:  numsp[s[i],i]=1, i=[0,nsver)
  GrB_Matrix numsp;
  GrB_Matrix_new(&numsp, GrB_INT32, n, nsver);
  GrB_Matrix_build(numsp,s,i_nsver,ones,nsver,GrB_PLUS_INT32);
  free(i_nsver); free(ones);

  // frontier: Holds the current frontier where values are path counts.
  // Initialized to out vertices of each source node in s.
  GrB_Matrix frontier;
  GrB_Matrix_new(&frontier, GrB_INT32, n, nsver);
  GrB_extract(frontier,numsp,GrB_NULL,A,GrB_ALL,n,s,nsver,desc_tsr);

  // sigma: stores frontier information for each level of BFS phase.  The memory
  // for an entry in sigmas is only allocated within the do-while loop if needed
  GrB_Matrix *sigmas = (GrB_Matrix*)malloc(sizeof(GrB_Matrix)*n);   // n is an upper bound on diameter

  int32_t   d = 0;                                       // BFS level number
  GrB_Index nvals = 0;                                   // nvals == 0 when BFS phase is complete

  // --------------------- The BFS phase (forward sweep) ---------------------------
  do {
    // sigmas[d](:,s) = d^th level frontier from source vertex s
    GrB_Matrix_new(&(sigmas[d]), GrB_BOOL, n, nsver);

    GrB_apply(sigmas[d],GrB_NULL,GrB_NULL,
              GrB_IDENTITY_BOOL,frontier,GrB_NULL);    // sigmas[d](:,:) = (Boolean) frontier
    GrB_eWiseAdd(numsp,GrB_NULL,GrB_NULL,
                 Int32Add,numsp,frontier,GrB_NULL);    // numsp += frontier (accum path counts)
    GrB_mxm(frontier,numsp,GrB_NULL,
            Int32AddMul,A,frontier,desc_tsr);          // f<!numsp> = A' +.* f (update frontier)
    GrB_Matrix_nvals(&nvals,frontier);                 // number of nodes in frontier at this level
    d++;
  } while (nvals);

  GrB_Monoid FP32Add;                                  // Monoid <float,+,0.0>
  GrB_Monoid_new(&FP32Add,GrB_PLUS_FP32,0.0f);
  GrB_Monoid FP32Mul;                                  // Monoid <float,*,1.0>
  GrB_Monoid_new(&FP32Mul,GrB_TIMES_FP32,1.0f);
  GrB_Semiring FP32AddMul;                             // Semiring <float,float,float,+,*,0.0>
  GrB_Semiring_new(&FP32AddMul,FP32Add,GrB_TIMES_FP32);

  // nspinv: the inverse of the number of shortest paths for each node and starting vertex.  |\label{line:nspinv}|
  GrB_Matrix nspinv;
  GrB_Matrix_new(&nspinv,GrB_FP32,n,nsver);
  GrB_apply(nspinv,GrB_NULL,GrB_NULL,
            GrB_MINV_FP32,numsp,GrB_NULL);             // nspinv = 1./numsp

  // bcu: BC updates for each vertex for each starting vertex in s
  GrB_Matrix bcu;
  GrB_Matrix_new(&bcu,GrB_FP32,n,nsver);
  GrB_assign(bcu,GrB_NULL,GrB_NULL,
             1.0f,GrB_ALL,n, GrB_ALL,nsver,GrB_NULL);  // filled with 1 to avoid sparsity issues

  // Descriptor used in the tally phase
  GrB_Descriptor desc_r;
  GrB_Descriptor_new(&desc_r);
  GrB_Descriptor_set(desc_r,GrB_OUTP,GrB_REPLACE);     // clear output before result is stored

  GrB_Matrix w;                                        // temporary workspace matrix
  GrB_Matrix_new(&w,GrB_FP32,n,nsver);

  // -------------------- Tally phase (backward sweep) --------------------
  for (int i=d-1; i>0; i--)  {
    GrB_eWiseMult(w,sigmas[i],GrB_NULL,
                  FP32Mul,bcu,nspinv,desc_r);          // w<sigmas[i]>=(1 ./ nsp).*bcu

    // add contributions by successors and mask with that BFS level's frontier
    GrB_mxm(w,sigmas[i-1],GrB_NULL,
            FP32AddMul,A,w,desc_r);                    // w<sigmas[i-1]> = (A +.* w)
    GrB_eWiseMult(bcu,GrB_NULL,GrB_PLUS_FP32,
                  FP32Mul,w,numsp,GrB_NULL);           // bcu += w .* numsp
  }

  // subtract "nsver" from every entry in delta (account for 1 extra value per bcu element)
  GrB_assign(*delta,GrB_NULL,GrB_NULL,
             -(float)nsver,GrB_ALL,n,GrB_NULL);        // fill with -nsver
  GrB_reduce(*delta,GrB_NULL,GrB_PLUS_FP32,
             GrB_PLUS_FP32,bcu,GrB_NULL);              // add all updates to -nsver

  // Release resources
  for(int i=0; i<d; i++) {
    GrB_free(&(sigmas[i]));
  }
  free(sigmas);

  GrB_free(&frontier);     GrB_free(&numsp);
  GrB_free(&nspinv);       GrB_free(&bcu);       GrB_free(&w);
  GrB_free(&desc_tsr);     GrB_free(&desc_r);
  GrB_free(&Int32AddMul);  GrB_free(&Int32Add);
  GrB_free(&FP32AddMul);   GrB_free(&FP32Add);   GrB_free(&FP32Mul);

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
    float values[] = {1., 1., 1., 1., 1., 1.,
                      1., 1., 1., 1., 1., 1.};

    GrB_Matrix graph;
    GrB_Vector bc;

    // Initialize a GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix_new(&graph, GrB_FP32, NUM_NODES, NUM_NODES);
    GrB_Matrix_build(graph, row_indices, col_indices, (float*)values, NUM_EDGES,
                     GrB_LOR);

    pretty_print_matrix_BOOL(graph, "GRAPH");
    GrB_Index sources[] = {0,1,2,3,4,5,6};

    BC_update(&bc, graph, sources, NUM_NODES);  // Batch BC from all nodes
    pretty_print_vector_FP64(bc, "BC (single batch)");
    GrB_free(&bc);

    GrB_Vector total_bc;
    GrB_Vector_new(&total_bc, GrB_FP32, NUM_NODES);
    for (GrB_Index src = 0; src < NUM_NODES; ++src)
    {
        BC(&bc, graph, src);
        GrB_eWiseAdd(total_bc, GrB_NULL, GrB_NULL, GrB_PLUS_FP32, bc, total_bc, GrB_NULL);
        pretty_print_vector_FP64(bc, "BC for 1 source");
        GrB_free(&bc);
    }
    pretty_print_vector_FP64(total_bc, "BC (accumulated)");

    // Cleanup
    GrB_free(&graph);
    GrB_free(&total_bc);
    GrB_finalize();

    return GrB_SUCCESS;
}
