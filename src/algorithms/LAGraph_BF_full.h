//------------------------------------------------------------------------------
// LAGraph_BF_full.c: Bellman-Ford single-source shortest paths, returns tree
//------------------------------------------------------------------------------

#pragma once

#include "../../deps/GraphBLAS/Include/GraphBLAS.h"

GrB_Info LAGraph_BF_full
(
	GrB_Vector *pd_output,      //the pointer to the vector of distance
	GrB_Vector *ppi_output,     //the pointer to the vector of parent
	GrB_Vector *ph_output,       //the pointer to the vector of hops
	const GrB_Matrix A,         //matrix for the graph
	const GrB_Index s           //given index of the source
);
