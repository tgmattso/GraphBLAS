#
# This file is part of the GraphBLAS tutorial materials,
# Copyright (c) 2018 Carnegie Mellon University and Intel Corporation.
#
# THIS SOFTWARE IS PROVIDED "AS IS," WITH NO WARRANTIES WHATSOEVER. CARNEGIE
# MELLON UNIVERSITY AND INTEL CORPORATION EXPRESSLY DISCLAIMS TO THE FULLEST
# EXTENT PERMITTED BY LAW ALL EXPRESS, IMPLIED, AND STATUTORY WARRANTIES,
# INCLUDING, WITHOUT LIMITATION, THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF PROPRIETARY RIGHTS.
#
# Released under a BSD (SEI)-style license, please see LICENSE.txt for
# full terms.
#
# Authors: Scott McMillan, Timothy G. Mattson
# 

# Project root not used yet
PROJECT_ROOT	= some_path_to/HPEC18_tutorial/code

INCDIR      = ../../include
LIBDIR      = ../../lib

include ../make.def

SOURCES 	= bfs_level.c betweenness_centrality.c maximal_independent_set.c \
		triangle_count.c page_rank.c

##  hpec_ex2_neighbors.c hpec_ex3_level_bfs.c hpec_ex4_bfs_parents.c hpec_ex34_parent_bfs.c

OBJECTS		= $(SOURCES:.c=.o)
TARGETS		= $(SOURCES:.c=.exe)
DSYM   		= $(SOURCES:.c=.exe.dSYM)
HEADERS		= $(wildcard *.h)
PCHS		= $(HEADERS:=.gch)

.PHONY		: all clean

all		: $(TARGETS)

%.s 		: %.c tutorial_utils.h
		$(CC) -S $(CCFLAGS) $<

%.o  		: %.c tutorial_utils.h
		$(CC) -c $(CCFLAGS) $<

%.exe  		: %.o
		$(CC) $^ -o $@ $(LIB)

clean		:
		/bin/rm -rf a.out $(TARGETS) $(OBJECTS) $(PCHS) $(DSYM) *~


squeaky		: clean
		/bin/rm -f $(SOURCES:%.c=%.d)


-include $(SOURCES:%.c=%.d)
