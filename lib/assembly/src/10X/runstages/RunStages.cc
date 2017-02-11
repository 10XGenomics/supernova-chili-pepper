// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#include "10X/runstages/RunStages.h"
#include "feudal/SubsetMasterVec.h"
#include "10X/Super.h"
#include "10X/Rescue.h"

void StageClosures(HyperBasevectorX& hb, vec<int>& inv, ReadPathVec& paths,
          VecULongVec& paths_index, vec<Bool>& dup, vec<Bool>& bad,
          vec<vec<int>>& all_closures, digraphE<vec<int>>& D,
          vec<int>& dinv)
{
     STAGE(Closures);
     // Make closures and the graph from them.
     cout << Date( ) << ": start making closures" << endl;
     MakeClosures( hb, inv, paths, paths_index, dup, bad, all_closures, True );
     MEM(before_destroy_paths_index);
     Destroy(paths_index);
     MEM(before_closures_to_graph);
     ClosuresToGraph( hb, inv, all_closures, D, dinv, True );
     Validate( hb, inv, D, dinv );
     MEM(after_closures_to_graph);
}


void StageEBC( HyperBasevectorX& hb, vec<int>& inv, VecULongVec& paths_index,
               vec<int32_t>& bc, vec<vec<int>>& ebc, VecIntVec& ebcx )
{
     STAGE(EBC);
     cout << Date( ) << ": computing ebc" << endl;
     ebc.resize( hb.E( ) );
     #pragma omp parallel for schedule(dynamic, 1000)
          for ( int e = 0; e < hb.E( ); e++ )
               ebc[e] = GetBarcodes( e, inv, paths_index, bc );

     {
          ebcx.resize( ebc.size( ) );
          for ( int64_t i = 0; i < (int64_t) ebc.size( ); i++ )
               {    ebcx[i].resize( ebc[i].size( ) );
                    for ( int j = 0; j < (int) ebc[i].size( ); j++ )
                         ebcx[i][j] = ebc[i][j];    }
     }
}

void StageDups( HyperBasevectorX& hb, vecbasevector& bases,
          ObjectManager<VecPQVec>& quals_om, ReadPathVec& paths,
          vec<int32_t>& bc, vec<Bool>& dup, vec<Bool>& bad, double& interdup_rate)

{
     STAGE(Dups);
//     VecPQVec& quals = quals_om.load_mutable( );
     quals_om.unload();
     VirtualMasterVec<PQVec> vmv( quals_om.filename() );
//     MarkDups( bases, quals, paths, bc, dup, interdup_rate );
     MarkDups( bases, vmv, paths, bc, dup, interdup_rate );
//     MarkBads( hb, bases, quals, paths, bad );
     MarkBads( hb, bases, vmv, paths, bad );
//     quals_om.unload( );
}

void StageTrim( HyperBasevectorX& hb, HyperBasevector& hbv, vec<int>& inv,
          vecbasevector& bases, ObjectManager<VecPQVec>& quals_om,
          ReadPathVec& paths, VecULongVec& paths_index, vec<int32_t>& bc )
{
     STAGE(Trim);
     quals_om.unload();
     MEM(trim_unload_quals_om);
     VirtualMasterVec<PQVec> vmv( quals_om.filename() );
     vec<int> dels;
     vec<Bool> dup;
     double interdup;
     MarkDups( bases, vmv, paths, bc, dup, interdup );
     VecPQVec& quals = quals_om.load_mutable( );
     MowLawn( hb, inv, bases, quals, paths, paths_index, bc, dup, interdup, dels );
     quals_om.unload();
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     hbv = HyperBasevector(hb);
     hbv.DeleteEdgesParallel(dels);
     Cleanup( hbv, inv, paths );
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     hb = HyperBasevectorX(hbv);
     cout << Date( ) << ": inverting paths index, mem usage = "
          << MemUsageGBString( ) << endl;
     paths_index.clear( );
     invert( paths, paths_index, hbv.E( ) );

     // Let's kill some hanging ends.

     vec<int> dels2;
     for ( int e = 0; e < hb.E( ); e++ )
     {    int re = inv[e];
          if ( paths_index[e].size( ) + paths_index[re].size( ) > 0 )
               continue;
          if ( hb.Kmers(e) > 200 ) continue; // not thought out!!
          // Bool alt = False;
          int v = hb.ToLeft(e), w = hb.ToRight(e);
          if ( hb.To(v).size( ) == 0 && hb.From(v).size( ) == 1 )
          {    dels2.push_back( e, re );    }
          else if ( hb.From(w).size( ) == 0 && hb.To(w).size( ) == 1 )
          {    dels2.push_back( e, re );    }    }

     UniqueSort(dels2);
     cout << Date( ) << ": deleting " << dels2.size( ) << " hanging ends"
          << endl;
     hbv.DeleteEdgesParallel(dels2);
     cout << Date( ) << ": cleaning up" << endl;
     Cleanup( hbv, inv, paths );
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     hb = HyperBasevectorX(hbv);
     paths_index.clear( );
     invert( paths, paths_index, hbv.E( ) );

}


void StageExtension( HyperBasevector& hbv, vec<int>& inv,
          vecbasevector& bases, ObjectManager<VecPQVec>& quals_om,
          ReadPathVec& paths, Bool const BACK_EXTEND )
{
     STAGE(Extension);
//     VecPQVec& qualsx = quals_om.load_mutable( );
     quals_om.unload();
     VirtualMasterVec<PQVec> vquals(quals_om.filename());
     cout << "Stage Extension, quals loaded, now current mem = " << MemUsageGBString( ) << endl;
     cout << "paths size before ExtendPathsNew " << paths.SizeSum() << endl;
     ExtendPathsNew( hbv, inv, bases, vquals, paths, BACK_EXTEND );
     cout << "Stage Extension, paths extended, now current mem = " << MemUsageGBString( ) << endl;
     cout << "paths size after ExtendPathsNew " << paths.SizeSum() << endl;
     quals_om.unload();
     cout << "after quals unload now current mem = " << MemUsageGBString( ) << endl;
}

void StagePatch(String const& dir, int const K, vecbasevector& bases, 
          ObjectManager<VecPQVec>& quals_om, HyperBasevector& hbv, 
          HyperBasevectorX& hb, ReadPathVec& paths, VecULongVec& paths_index,
          vec<int>& inv, const vec<Bool>& dup, vec<Bool>& bad, 
          vec<DataSet>& datasets, vec<int32_t>& bc, const int max_width, 
          Bool ONE_GOOD, vec<basevector>& closures, vec<pair<int,int>>& pairs, 
          Bool const CG2, const Bool STACKSTER,
          const Bool STACKSTER_ALT, const Bool RESCUE )
{
     // Get started.

     STAGE(Patch);
     ForceAssertEq( dup.size()*2, bases.size() );
     MEM(patch_start);
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     hb = HyperBasevectorX(hbv);
     MEM(hbv_copy);

     // Rescue kmers.

     if (RESCUE) Rescue( 0, bases, bc, datasets, hb, inv, paths, closures );

     // Proceed with the rest of the patching.

     VirtualMasterVec<PQVec> vmv( quals_om.filename() );
     quals_om.unload();
//     VecPQVec& quals = quals_om.load_mutable( );
//     MarkBads( hb, bases, quals, paths, bad );
     MarkBads( hb, bases, vmv, paths, bad );

     // Look for edge pairs that we should try to close.

     cout << Date( ) << ": finding edge pairs" << endl;
     FindEdgePairs(
          hb, inv, paths, paths_index, bad, pairs, datasets, bc, ONE_GOOD );
     MEM(edge_pairs);

     // Determine which read ids will be used

     cout << Date() << ": determining subset of interest" << endl;
    vec<uint64_t> read_ids;
#if 1
     for ( int pi = 0; pi < pairs.isize( ); pi++ )
     {    int e1 = pairs[pi].first, e2 = pairs[pi].second;
          vec<int> es = { e1, inv[e2] };
          for ( int ei = 0; ei < 2; ei++ )
          {    int e = es[ei];
               for ( int pass = 1; pass <= 2; pass++ )
               {    const int f = ( pass == 1 ? e : inv[e] );
                    for ( int i = 0; i < (int) paths_index[f].size( ); i++ )
                    {    
                         uint64_t id = paths_index[f][i];
                         read_ids.push_back(id);
                         uint64_t id2 = ( id % 2 == 0 ? id+1 : id-1 );
                         read_ids.push_back(id2);
                    }
               }
          }
     }
     cout << Date() << ": uniquesort" << endl;
     UniqueSort(read_ids);
#else
     for (size_t i = 0; i < paths.size(); ++i )   // make this a NO-OP for testing
          read_ids.push_back(i);
#endif
     cout << Date() << ": subset is size " << read_ids.size() << endl;

     cout << Date() << ": Destroying paths and populating subsets from disk" << endl;
     Destroy(paths);
     MEM(destroy_paths);

     SubsetMasterVec<PQVec> quals_subset( quals_om.filename(), read_ids );
     SubsetMasterVec<ReadPath> paths_subset( dir+"/a.paths", read_ids );

     MEM(after_populate_subset);

#if 0
     // TEMP - check that the NO-OP worked
     cout << Date() << ": checking paths and quals" << endl;
     ForceAssertEq(quals_subset.size(), quals.size() );
     ForceAssertEq(paths_subset.size(), paths.size() );
     for ( size_t i = 0; i < quals_subset.size(); ++i ) {
          if ( paths_subset[i] !=  paths[i] ) FatalErr( "paths bad" );
          for ( size_t j = 0; j < quals_subset[i].size(); ++j )  {
               qualvector q1, q2;
               quals_subset[i].unpack(&q1);
               quals[i].unpack(&q2);
               if ( q1 != q2 ) FatalErr("quals bad" );
          }
     }
     cout << Date() << ": done checking paths and quals" << endl;
#endif


     // Traverse pairs.

     cout << Date( ) << ": start traversing " << ToStringAddCommas( pairs.size( ) )
          << " pairs" << endl;
     double pclock = WallClockTime( );
     const int batches = 1000;
     vec< vec<basevector> > closuresi(batches);
     int64_t NP = pairs.size( );
     vec<int> kmers( hb.E( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          kmers[e] = hb.Kmers(e);
//     auto& quals = quals_om.load_mutable();       // TO-DO: REMOVE ME 
     #pragma omp parallel for
     for ( int bi = 0; bi < batches; bi++ )
     {    for ( int pi = bi * NP / batches; pi < (bi+1) * NP / batches; pi++ )
          {    int e1 = pairs[pi].first, e2 = pairs[pi].second;
               if ( STACKSTER )
               {    vec<basevector> edges = { hb.O(e1), hb.O(e2) };
                    const int VERBOSITY = 0;
                    vec<basevector> f;
                    const Bool EXP = False;
                    vec<int> trim;
                    Stackster( e1, e2, edges, bases, quals_subset, K, datasets, 
                         kmers, inv, dup, paths_subset, paths_index, f, trim,
                         VERBOSITY, STACKSTER_ALT, EXP );
                    closuresi[bi].append(f);    }
               Bool verbose = False;
               vec<basevector> f;
               if ( CG2 ) CloseGap2( hb, inv, bases, quals_subset, paths_subset, paths_index,
                    e1, e2, f, verbose, pi, max_width );
               else CloseGap( hb, inv, bases, quals_subset, paths_subset, paths_index,
                    e1, e2, f, verbose, pi, max_width );
               {    closuresi[bi].append(f);    }    }    }

     MEM(traverse_pairs);
     Destroy(paths_index);
     MEM(destroy_paths_index);
     quals_om.unload( );
     MEM(quals_unload);
     Destroy(paths_subset);
     MEM(destroyed_paths_subset);
     Destroy(quals_subset);
     MEM(destroyed_quals_subset);
     Destroy(bases);
     MEM(destroyed_bases);

     for ( int bi = 0; bi < batches; bi++ )
          closures.append( closuresi[bi] );
     cout << Date( ) << ": found " << closures.size( ) << " closures" << endl;
     cout << TimeSince(pclock) << " used closing pairs" << endl;
     // cout << Date( ) << ": sorting closures" << endl;
     // ParallelSort(closures);

     // Insert gap patches into assembly.

     ForceAssertEq( K, hbv.K( ) );
     HyperBasevector hb3;
     vec<vec<int>> to3( hbv.E( ) );
     vec<int> left3( hbv.E( ) );

     MEM(start_build);
     {    ReadPathVec allx_paths;
          vecbasevector allx;
          BuildAll( allx, hbv, closures.size( ) );
          MEM(build_all);
          allx.append( closures.begin( ), closures.end( ) );
          cout << Date( ) << ": building hb2" << endl;
          cout << "memory in use now = "
               << ToStringAddCommas( MemUsageBytes( ) ) << endl;
          double clock2 = WallClockTime( );
          const int coverage = 4;
          MEM(start_build_big);
          buildBigKHBVFromReads( K, allx, coverage, &hb3, &allx_paths );
          MEM(build_big);
          cout << Date( ) << ": back from buildBigKHBVFromReads" << endl;

          // build to3 and left3 from allx_paths

          for ( int i = 0; i < hb.E( ); ++i )
          {    for ( auto const& p : allx_paths[i] )
                    to3[i].push_back( p );
               left3[i] = allx_paths[i].getFirstSkip();    }

          cout << TimeSince(clock2) << " used in new stuff 2 test" << endl;    }

     cout << "peak mem usage = " << PeakMemUsageGBString( ) << endl;
     MEM(before_load_paths);
     paths.ReadAll( dir + "/a.paths" );
     MEM(load_paths);
     double clock3 = WallClockTime( );
     TranslatePaths( paths, hb3, to3, left3 );
     Destroy(to3), Destroy(left3);
     hbv = hb3;
     hbv.Involution(inv);
     Validate( hbv, inv, paths );
}




void StageBuildGraph( int const K, vecbasevector& bases, ObjectManager<VecPQVec>& quals_om,
          int MIN_QUAL, int MIN_FREQ, int MIN_BC, vec<int32_t> const& bc, int64_t bc_start,
          std::string const GRAPH, double const GRAPHMEM,
          String work_dir, HyperBasevector& hbv, ReadPathVec& paths, VecULongVec& paths_index, vec<int>& inv)
{
     STAGE(BuildGraph);

     cout << Date() << ": barcoded *datatypes* start at " << bc_start << endl;

     MEM(before_graph_creation);
     if ( K != 48 ) FatalErr("remember that we have changes for K=48 not propagated to other K values.");
     buildReadQGraph48(work_dir, GRAPH, bases, quals_om, False, False, MIN_QUAL, MIN_FREQ, bc_start, MIN_BC, &bc,
                              .75, 0, "", True, False, &hbv, &paths, GRAPHMEM, False );
     cout << Date( ) << ": back from buildReadQGraph" << endl;
     cout << Date( ) << ": memory in use = " << MemUsageGBString( )
          << ", peak = " << PeakMemUsageGBString( ) << endl;
     hbv.Involution(inv);

     // Report assembly statistic.

     if ( hbv.E( ) == 0 )
     {    cout << "\nAssembly is empty, giving up.\n" << endl;
               Scram(1);    }
     vec<int> len( hbv.E( ) );
     for ( int e = 0; e < hbv.E( ); e++ )
          len[e] = hbv.Kmers(e);
     Sort(len);
     cout << "N50 edge length = " << N50(len) << endl;

     // Write files.

#if 0          // should be able to retire FixPaths now
     auto bads = FixPaths( hbv, paths ); // needed? - let's find out... see below
     if ( bads ) cout << Date() << ": FixPaths truncated " << bads << " BAD paths" << endl;
#endif

     cout << Date( ) << ": inverting paths index, mem usage = "
          << MemUsageGBString( ) << endl;
     invert( paths, paths_index, hbv.E( ) );
}


