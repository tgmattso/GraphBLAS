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

SOURCES 	= BuildAdjMatTuple.c \
		BuildAdjMatIndVec.c \
		matvec.c \
		matvecTrans.c \
		matvecTransIter.c \
		matvecTransIterVisited.c \
		matvecTransIterVisitedExitFlag.c \
		matvecTransIterLevels.c \
		bfs_level.c \
		connected_components.c \
		AnalyzeGraph_final.c

OBJECTS		= $(SOURCES:.c=.o)
TARGETS		= $(SOURCES:.c=.exe)
DSYM   		= $(SOURCES:.c=.exe.dSYM)
HEADERS		= $(wildcard $(INCDIR)/*.h)
PCHS		= $(HEADERS:=.gch)

.PHONY		: all clean

all		: $(TARGETS)

%.s 		: %.c $(HEADERS)
		$(CC) -S $(CCFLAGS) $<

%.o  		: %.c $(HEADERS)
		$(CC) -c $(CCFLAGS) $<

%.exe  		: %.o
		$(CC) $^ -o $@ $(LIB)

clean		:
		/bin/rm -rf a.out $(TARGETS) $(OBJECTS) $(PCHS) $(DSYM) *~


squeaky		: clean
		/bin/rm -f $(SOURCES:%.c=%.d)

.SECONDARY	: $(OBJECTS)

-include $(SOURCES:%.c=%.d)
