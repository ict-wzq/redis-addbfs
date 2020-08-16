#pragma once

#include "../../deps/GraphBLAS/Include/GraphBLAS.h"

GrB_Info bfs6               // BFS of a graph (using unary operator)
(
    GrB_Vector *v_output,   // v [i] is the BFS level of node i in the graph
    const GrB_Matrix A,     // input graph, treated as if boolean in semiring
    GrB_Index s             // starting node of the BFS
);
