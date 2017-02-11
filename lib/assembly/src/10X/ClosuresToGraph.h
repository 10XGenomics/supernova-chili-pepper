// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_CLOSURES_TO_GRAPH_H
#define TENX_CLOSURES_TO_GRAPH_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"

// Convert a digraphE<int> into a digraphE<vec<int>> in which unneeded vertices
// have been removed.

void Vectorify(

     // input: pointer to digraphE<int> -- GETS DELETED MIDSTREAM!!

     digraphE<int>* D0p,

     // input and output:

     vec<int>& dinv,

     // output:

     digraphE<vec<int>>& D,

     // control:

     const Bool verbose, const Bool single, const Bool use_inv );

void ClosuresToGraph( const HyperBasevectorX& hb, const vec<int>& inv, 
     vec<vec<int>>& all_closures, digraphE<vec<int>>& D, vec<int>& dinv,
     const Bool verbose );

#endif
