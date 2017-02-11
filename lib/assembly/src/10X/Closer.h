// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_CLOSER_H
#define TENX_CLOSER_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

template< class VP, class VPI > void Closer( const HyperBasevectorX& hb, 
     const vec<int>& inv, VP& paths, VPI& paths_index, const vec<Bool>& dup, 
     const vec<Bool>& bad, vec<vec<int>>& all_closures,
     int64_t& opcount, const Bool GLOBAL, vec<int64_t>& pids, const int verbosity );

#endif
