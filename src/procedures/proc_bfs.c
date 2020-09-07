/*
* Copyright 2018-2020 Redis Labs Ltd. and Contributors
*
* This file is available under the Redis Labs Source Available License Agreement
*/

#include "proc_bfs.h"
#include "proc_ctx.h"
#include "../query_ctx.h"
#include "../datatypes/array.h"
#include "../graph/entities/node.h"
#include "../graph/graphcontext.h"
#include "../algorithms/bfs6.h"

// CALL algo.Bfs(startNode, label, relationship)

typedef struct {
	int n;                          // Number of nodes to do Bfs.
	Graph *g;                       // Graph.
	Node *startNode;                // Start Node to do Bfs.
	GrB_Index *mappings;            // Mappings between extracted matrix rows and node ids.
	GrB_Matrix M;					// relation Matrix whose nodes are labeled 
	SIValue *output;                // Array with 4 entries ["node", node, "level", level].
	bool depleted;					// True if BFS has already been performed for this node.
} BfsContext;

ProcedureResult Proc_BfsInvoke(ProcedureCtx *ctx, const SIValue *args) {
	if(array_len((SIValue *)args) != 3) return PROCEDURE_ERR;
	if(SI_TYPE(args[0]) != T_NODE || SI_TYPE(args[1]) != T_STRING || SI_TYPE(args[2]) != T_STRING) 
		return PROCEDURE_ERR;

	Node *startNode = (Node *)(args[0].ptrval);	
	const char *label = args[1].stringval;
	const char *relation = args[2].stringval;

	Graph *g = QueryCtx_GetGraph();
	GrB_Index n = 0;
	GrB_Index *mappings = NULL; // Mappings, array for returning row indices of tuples.
	GrB_Matrix reduced = GrB_NULL;

	// Setup context.
	BfsContext *pdata = rm_malloc(sizeof(BfsContext) * 1);
	pdata->n = n;
	pdata->g = g;
	pdata->startNode = startNode;
	pdata->mappings = mappings;
	pdata->M = reduced;
	pdata->output = array_new(SIValue, 4);
	pdata->output = array_append(pdata->output, SI_ConstStringVal("node_id"));
	pdata->output = array_append(pdata->output, SI_NullVal()); // Place holder.
	pdata->output = array_append(pdata->output, SI_ConstStringVal("level"));
	pdata->output = array_append(pdata->output, SI_NullVal()); // Place holder.
	pdata->depleted = false;
	ctx->privateData = pdata;

	// Get label matrix.
	Schema *s = NULL;
	GrB_Matrix l = NULL;
	GrB_Matrix r = NULL;
	GraphContext *gc = QueryCtx_GetGraphCtx();
	s = GraphContext_GetSchema(gc, label, SCHEMA_NODE);
	if(!s) return PROCEDURE_OK;// Failed to find schema, first step will return NULL.
	l = Graph_GetLabelMatrix(g, s->id);

	// Get relation matrix.
	s = GraphContext_GetSchema(gc, relation, SCHEMA_EDGE);
	if(!s) return PROCEDURE_OK;	
	r = Graph_GetRelationMatrix(g, s->id);

	//Get matrix for bfs
	GrB_Index rows = Graph_RequiredMatrixDim(g);
	GrB_Index cols = rows;
	printf("rows = %ld\n", rows);
	assert(GrB_Matrix_nvals(&n, l) == GrB_SUCCESS);
	assert(GrB_Matrix_new(&reduced, GrB_BOOL, n, n) == GrB_SUCCESS);
	printf("n = %ld\n", n);
	if(n != rows) {
		mappings = rm_malloc(sizeof(GrB_Index) * n);
		assert(GrB_Matrix_extractTuples_BOOL(mappings, GrB_NULL, GrB_NULL, &n, l) == GrB_SUCCESS);
		assert(GrB_extract(reduced, GrB_NULL, GrB_NULL, r, mappings, n, mappings, n,
						   GrB_NULL) == GrB_SUCCESS);//extract submatrix from r
		pdata->M = reduced;
	} else {
		pdata->M = r;
	}

	// Update context.
	pdata->n = n;
	pdata->mappings = mappings;

	return PROCEDURE_OK;
}

SIValue *Proc_BfsStep(ProcedureCtx *ctx) {
	assert(ctx->privateData);
	BfsContext *pdata = (BfsContext *)ctx->privateData;
	
	// Return NULL if this source has already been mapped or there are no connected nodes.
	if(pdata->depleted || pdata->n == 0) return NULL;

	Node *s = pdata->startNode;
	NodeID s_id = ENTITY_GET_ID(s);
	printf("s_id=%ld\n",s_id);
	GrB_Vector output = GrB_NULL;    // Pointer to the vector of level

	assert(bfs6(&output, pdata->M, s_id) == GrB_SUCCESS);

	//Get number of entries in bfs output
	GrB_Index nvals;
	GrB_Vector_nvals(&nvals, output);
	SIValue nodes = SI_Array(nvals);
	SIValue levels = SI_Array(nvals);
	printf("nvals of output = %ld\n", nvals);

	//set result
	GrB_Index *node_ids = rm_malloc(nvals * sizeof(GrB_Index));
	int32_t *level = rm_malloc(nvals * sizeof(int32_t));
	assert(GrB_Vector_extractTuples_INT32(node_ids, level, &nvals, output) == GrB_SUCCESS);
	for (GrB_Index i = 0; i < nvals; i++) {
		NodeID node_id = (pdata->mappings) ? pdata->mappings[i] : i;
		SIArray_Append(&nodes, SI_LongVal(node_id));
		SIArray_Append(&levels, SI_LongVal(level[i]));
		printf("%ld : node = %ld, level = %d\n", i, node_id, level[i]);
	}	
	pdata->output[1] = nodes;
	pdata->output[3] = levels;

	pdata->depleted = true; // Mark that this node has been mapped.
	return pdata->output;
}

ProcedureResult Proc_BfsFree(ProcedureCtx *ctx) {
	// Clean up.
	if(ctx->privateData) {
		BfsContext *pdata = ctx->privateData;
		GrB_free(&pdata->M);
		array_free(pdata->output);
		rm_free(pdata->mappings);
		rm_free(ctx->privateData);
	}

	return PROCEDURE_OK;
}

ProcedureCtx *Proc_BfsCtx() {
	void *privateData = NULL;

	ProcedureOutput **outputs = array_new(ProcedureOutput *, 2);
	ProcedureOutput *output_node = rm_malloc(sizeof(ProcedureOutput));
	ProcedureOutput *output_level = rm_malloc(sizeof(ProcedureOutput));

	output_node->name = "node_id";
	output_node->type = T_ARRAY;

	output_level->name = "level";
	output_level->type = T_ARRAY;

	outputs = array_append(outputs, output_node);
	outputs = array_append(outputs, output_level);
	ProcedureCtx *ctx = ProcCtxNew("algo.BFS",
								   3,
								   outputs,
								   Proc_BfsStep,
								   Proc_BfsInvoke,
								   Proc_BfsFree,
								   privateData,
								   true);
	return ctx;
}
