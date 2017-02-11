// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_N50P_H
#define TENX_N50P_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"

int64_t N50P( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<int>& dinv, const vecbasevector& G, 
     const vec<vec<pair<int,int>>>& locs, const vec<Bool>& keep, 
     String& report, String& details, double& errw );

#endif
