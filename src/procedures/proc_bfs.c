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
} BfsContext;

ProcedureResult Proc_BfsInvoke(ProcedureCtx *ctx, const SIValue *args) {
	if(array_len((SIValue *)args) != 3) return PROCEDURE_ERR;
	if(!(SI_TYPE(args[1]) & SI_TYPE(args[2]) & T_STRING)) return PROCEDURE_ERR;

	Node *startNode = (Node *)(args[0].ptrval);	
	const char *label = args[1].stringval;
	const char *relation = args[2].stringval;

	GrB_Index n = 0;
	Schema *s = NULL;
	GrB_Matrix l = NULL;
	GrB_Matrix r = NULL;
	GrB_Index *mappings = NULL; // Mappings, array for returning row indices of tuples.
	GrB_Matrix reduced = GrB_NULL;
	Graph *g = QueryCtx_GetGraph();
	GraphContext *gc = QueryCtx_GetGraphCtx();

	// Setup context.
	BfsContext *pdata = rm_malloc(sizeof(BfsContext));
	pdata->g = g;

	// Get label matrix.
	s = GraphContext_GetSchema(gc, label, SCHEMA_NODE);
	if(!s) return PROCEDURE_OK;
	l = Graph_GetLabelMatrix(g, s->id);

	// Get relation matrix.
	s = GraphContext_GetSchema(gc, relation, SCHEMA_EDGE);
	if(!s) return PROCEDURE_OK;
	r = Graph_GetRelationMatrix(g, s->id);

	//Get matrix for bfs
	GrB_Index rows = Graph_RequiredMatrixDim(g);
	GrB_Index cols = rows;
	assert(GrB_Matrix_nvals(&n, l) == GrB_SUCCESS);
	assert(GrB_Matrix_new(&reduced, GrB_BOOL, n, n) == GrB_SUCCESS);
	
	if(n != rows) {
		mappings = rm_malloc(sizeof(GrB_Index) * n);
		assert(GrB_Matrix_extractTuples_BOOL(mappings, GrB_NULL, GrB_NULL, &n, l) == GrB_SUCCESS);
		assert(GrB_extract(reduced, GrB_NULL, GrB_NULL, r, mappings, n, mappings, n,
						   GrB_NULL) == GrB_SUCCESS);
	} else {
		/* There no need to perform extraction as `r` dimension NxN
		 * is the same as the number of entries in `l` which means
		 * all connections described in `r` connect nodes of type `l`
		 * Unfortunately we still need to type cast `r` into a boolean matrix. */
		GrB_Descriptor desc;
		GrB_Descriptor_new(&desc);
		GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
		assert(GrB_transpose(reduced, GrB_NULL, GrB_NULL, r, desc) == GrB_SUCCESS);
		GrB_free(&desc);
	}

	// Update context.
	pdata->n = n;
	pdata->M = reduced;
	pdata->mappings = mappings;

	pdata->output = array_new(SIValue, 4);
	pdata->output = array_append(pdata->output, SI_ConstStringVal("node_id"));
	pdata->output = array_append(pdata->output, SI_NullVal()); // Place holder.
	pdata->output = array_append(pdata->output, SI_ConstStringVal("level"));
	pdata->output = array_append(pdata->output, SI_NullVal()); // Place holder.

	ctx->privateData = pdata;
	return PROCEDURE_OK;
}

SIValue *Proc_BfsStep(ProcedureCtx *ctx) {
	assert(ctx->privateData);

	BfsContext *pdata = (BfsContext *)ctx->privateData;
	int32_t v = 0;
	GrB_Index n = pdata->n;
	SIValue node = SIArray_New(n);
	SIValue level = SIArray_New(n);
	Node *s = pdata->startNode;
	NodeID s_id = ENTITY_GET_ID(s);
	GrB_Vector output = GrB_NULL;    // Pointer to the vector of level

	GrB_Info info = bfs6(&output, pdata->M, s_id);

	if(info != GrB_SUCCESS) {
		//Failed to run bfs6 , return NULL.
		return NULL;
	}

	//set result
	for (GrB_Index i = 0; i < n; i++) {
		NodeID node_id = (pdata->mappings) ? pdata->mappings[i] : i;
		SIArray_Append(&node, SI_LongVal(node_id));

		info = GrB_Vector_extractElement_INT32(&v, output, i);
		if(info != GrB_SUCCESS) {
			//Failed to extract element from output , return NULL.
			return NULL;
		}
		SIArray_Append(&level, SI_LongVal(c));
	}

	//Graph_GetNode(pdata->g, node_id, &pdata->node);
	pdata->output[1] = node;
	pdata->output[3] = level;

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

ProcedureCtx *Proc_PagerankCtx() {
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
	ProcedureCtx *ctx = ProcCtxNew("algo.Bfs",
								   3,
								   outputs,
								   Proc_BfsStep,
								   Proc_BfsInvoke,
								   Proc_BfsFree,
								   privateData,
								   true);
	return ctx;
}