// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef _10X_RUNSTAGES_H
#define _10X_RUNSTAGES_H
#include "10X/DfTools.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "Intvector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "graph/DigraphTemplate.h"
#include "kmers/BigKPather.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/UnibaseUtils.h"
#include "paths/long/BuildReadQGraph40.h"
#include "paths/long/BuildReadQGraph48.h"
#include "paths/long/BuildReadQGraph60.h"
#include "paths/long/ReadPath.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/large/Repath.h"
#include "system/HostName.h"
#include "10X/Closomatic.h"
#include "10X/ClosuresToGraph.h"
#include "10X/DfTools.h"
#include "10X/Extend.h"
#include "10X/Heuristics.h"
#include "10X/Lawnmower.h"
#include "10X/SecretOps.h"
#include "10X/Stackster.h"
#include "10X/WriteFiles.h"
#include "10X/MaybeWriter.h"
#include "10X/astats/RefLookup.h"
#include "10X/Predup.h"


void StageClosures(HyperBasevectorX& hb, vec<int>& inv, ReadPathVec& paths,
          VecULongVec& paths_index, vec<Bool>& dup, vec<Bool>& bad,
          vec<vec<int>>& all_closures, digraphE<vec<int>>& D,
          vec<int>& dinv);

void StageEBC( HyperBasevectorX& hb, vec<int>& inv, VecULongVec& paths_index,
               vec<int32_t>& bc, vec<vec<int>>& ebc, VecIntVec& ebcx );

void StageDups( HyperBasevectorX& hb, vecbasevector& bases,
          ObjectManager<VecPQVec>& quals_om, ReadPathVec& paths,
          vec<int32_t>& bc, vec<Bool>& dup, vec<Bool>& bad, double& interdup_rate);

void StageTrim( HyperBasevectorX& hb, HyperBasevector& hbv, vec<int>& inv,
          vecbasevector& bases, ObjectManager<VecPQVec>& quals_om,
          ReadPathVec& paths, VecULongVec& paths_index, vec<int32_t>& bc );

void StageExtension( HyperBasevector& hbv, vec<int>& inv,
          vecbasevector& bases, ObjectManager<VecPQVec>& quals_om,
          ReadPathVec& paths, Bool const BACK_EXTEND );

void StagePatch(String const& dir, int const K, vecbasevector& bases, ObjectManager<VecPQVec>& quals_om,
          HyperBasevector& hbv, HyperBasevectorX& hb, ReadPathVec& paths,
          VecULongVec& paths_index, vec<int>& inv, const vec<Bool>& dup,
          vec<Bool>& bad, vec<DataSet>& datasets,
          vec<int32_t>& bc, const int max_width, Bool ONE_GOOD, vec<basevector>& closures,
          vec<pair<int,int>>& pairs, Bool const CG2, const Bool STACKSTER,
          const Bool STACKSTER_ALT, const Bool RESCUE );

void StageBuildGraph( int const K, vecbasevector& bases, ObjectManager<VecPQVec>& quals_om,
          int MIN_QUAL, int MIN_FREQ, int MIN_BC, vec<int32_t> const& bc, int64_t bc_start,
          std::string const GRAPH, double const GRAPHMEM,
          String work_dir, HyperBasevector& hbv, ReadPathVec& paths, VecULongVec& paths_index, vec<int>& inv);

void StageEBC( HyperBasevectorX& hb, vec<int>& inv, VecULongVec& paths_index,
               vec<int32_t>& bc, vec<vec<int>>& ebc, VecIntVec& ebcx );

struct RunStages {
     RunStages( String const& START, vector<const char*> const& init ) {
          size_t i=0;
          for ( auto const& stage : init ) _stages[string(stage)]=++i;
          cout << "START=" << START << endl;
          try {
               _start = _stages.at(START);
          } catch ( std::out_of_range& )  {
               ostringstream s;
               s << "unknown START stage: " << START;
               FatalErr( s.str() );
          }
     }
     bool at_or_after( string const& here ) {
          try { return _stages.at( here ) <= _start; }
          catch ( std::out_of_range& ) { FatalErr("bad stage in at_or_before()"); }
     }
     bool at_or_before( string const& here ) {
          try { return _stages.at( here ) >= _start; }
          catch ( std::out_of_range& ) { FatalErr("bad stage in at_or_before()"); }
     }
     bool at( string const& here ) {
          try { return _stages.at( here ) == _start; }
          catch ( std::out_of_range& ) { FatalErr("bad stage in at()" ); }
     }
private:
     size_t _start;
     map< string, size_t > _stages;
};

#define STAGE(x) { cout << Date() << ": ** start stage " << #x << ", mem = " << MemUsageGBString() << ", peak = " << PeakMemUsageGBString() << endl; }




#endif /* _10X_RUNSTAGES_H */
