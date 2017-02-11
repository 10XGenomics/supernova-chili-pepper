// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/ClosuresToGraph.h"
#include "10X/WDG.h"
#include "10X/Super.h"

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

     const Bool verbose, const Bool single, const Bool use_inv )
{    
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     vec<int> to_left, to_right, edges0;
     vec<int> dinv2;
     {    digraphE<int>& D0 = *D0p;
          if (verbose)
          {    cout << Date( ) << ": converting, peak mem = "
                    << PeakMemUsageGBString( ) << endl;    }
          D0.ToLeft(to_left), D0.ToRight(to_right);
          vec<Bool> used( D0.E( ), False );
          for ( int e = 0; e < D0.E( ); e++ )
          {    Bool lverbose = False;
               if ( used[e] ) continue;
               if ( lverbose ) cout << "past continue for e=" << e << endl;
               vec<int> x = {e};
               // used[e] = True; // true but pointless
               Bool circle = False;
               while(1)
               {    int e = x.back( );
                    int v = to_right[e];
                    if ( !D0.From(v).solo( ) || !D0.To(v).solo( ) ) break;
                    int f = D0.IFrom( v, 0 );
                    if ( f == x[0] )
                    {    circle = True;
                         break;    }
                    x.push_back(f);
                    used[f] = True;    }
               if ( lverbose ) 
               {    cout << "grew forward for e=" << e << ", x=" << printSeq(x) 
                         << endl;
                    cout << "circle is " << (circle?"True":"False") << endl;    }
               if ( !circle )
               {    x.ReverseMe( );
                    while(1)
                    {    int e = x.back( );
                         int v = to_left[e];
                         if ( !D0.From(v).solo( ) || !D0.To(v).solo( ) ) break;
                         int f = D0.ITo( v, 0 );
                         if ( f == x[0] )
                         {    circle = True;
                              break;    }
                         x.push_back(f);
                         used[f] = True;    }
                    x.ReverseMe( );    }
               if ( lverbose )
               {    cout << "grew backward for e=" << e << ", x=" << printSeq(x) 
                         << endl;
                    cout << "circle is " << (circle?"True":"False") << endl;    }
               int E = D.E( );
               D.EdgesMutable( ).push_back(x);
               if (use_inv)
               {    vec<int> rx( x.size( ) );
                    for ( int i = 0; i < x.isize( ); i++ )
                    {    rx[ x.isize( ) - i - 1 ] = dinv[ x[i] ];
                         if ( lverbose ) PRINT3( i, x[i], dinv[x[i]] );    }
                    if ( lverbose ) cout << "rx is " << printSeq(rx) << endl;
                    for ( int i = 0; i < rx.isize( ); i++ )
                         used[ rx[i] ] = True;
                    int idx = E;    // index of the one we pushed on
                    if ( x != rx ) 
                    {    D.EdgesMutable().push_back(rx);
                         dinv2.push_back( idx+1 );    }
                    dinv2.push_back( idx );
                    AssertEq(D.E(), dinv2.isize() );    }    }
          if (verbose)
          {    cout << Date( ) << ": copying edges, peak mem = "
                    << PeakMemUsageGBString( ) << endl;    }
          // Note that the following copy seems "theoretically" unnecessary:
          edges0 = D0.Edges( );    
          delete D0p;    }
     if (use_inv) dinv = dinv2;
     if (verbose)
     {    cout << Date( ) << ": making verts, peak mem = "
               << PeakMemUsageGBString( ) << endl;    }
     vec<int> verts( 2 * D.E( ) );
     #pragma omp parallel for schedule(dynamic, 10000) num_threads(nthreads)
     for ( int i = 0; i < D.E( ); i++ )
     {    verts[2*i] = to_left[ D.O(i).front( ) ];
          verts[2*i+1] = to_right[ D.O(i).back( ) ];    }
     if (verbose)
     {    cout << Date( ) << ": sorting verts, peak mem = "
               << PeakMemUsageGBString( ) << endl;    }
     if (single) UniqueSort(verts);
     else ParallelUniqueSort(verts);
     int NV = verts.size( );
     if (verbose)
     {    cout << Date( ) << ": resizing, peak mem = " << PeakMemUsageGBString( )
               << endl;    }
     D.FromMutable( ).resize(NV), D.ToMutable( ).resize(NV);
     D.FromEdgeObjMutable( ).resize(NV), D.ToEdgeObjMutable( ).resize(NV);
     if (verbose)
     {    cout << Date( ) << ": filling out graph, peak mem = "
               << PeakMemUsageGBString( ) << endl;    }
     for ( int e = 0; e < D.E( ); e++ )
     {    int v = BinPosition( verts, to_left[ D.O(e).front( ) ] );
          int w = BinPosition( verts, to_right[ D.O(e).back( ) ] );
          D.FromMutable(v).push_back(w), D.ToMutable(w).push_back(v);
          D.FromEdgeObjMutable(v).push_back(e);
          D.ToEdgeObjMutable(w).push_back(e);    }
     Destroy(to_left), Destroy(to_right), Destroy(verts);
     if (verbose)
     {    cout << Date( ) << ": fixing edges, peak mem = "
               << PeakMemUsageGBString( ) << endl;    }
     #pragma omp parallel for schedule(dynamic, 10000) num_threads(nthreads)
     for ( int e = 0; e < D.E( ); e++ )
     {    for ( int j = 0; j < D.O(e).isize( ); j++ )
               D.OMutable(e)[j] = edges0[ D.O(e)[j] ];    }
     Destroy(edges0);
     if (verbose)
     {    cout << Date( ) << ": completing D construction, peak mem = "
               << PeakMemUsageGBString( ) << endl;    }
     #pragma omp parallel for schedule(dynamic, 10000) num_threads(nthreads)
     for ( int v = 0; v < D.N( ); v++ )
     {    SortSync( D.FromMutable(v), D.FromEdgeObjMutable(v) );
          SortSync( D.ToMutable(v), D.ToEdgeObjMutable(v) );    }

     // Done.

     if (verbose)
          cout << Date( ) << ": final graph has " << D.E( ) << " edges" << endl;    }

template <class V>
bool ValidatePath( V p, HyperBasevector const& hb,
          vec<int> const& to_left, vec<int> const& to_right )
{
     for ( int i = 1; i < p.size(); ++i )  {
          if ( to_left[p[i]] != to_right[p[i-1]] ) return false;
     }

     return true;
}

void ValidateD0( HyperBasevector const& hb, vec<int> const& inv, digraphE<int> const& D0, vec<int> const& dinv )
{
     // validate dinv
     cout << Date() << ": checkout dinv consistency with inv" << endl;
     ForceAssertEq( dinv.size(), D0.Edges().size() );
     for ( int i = 0; i < dinv.isize(); ++i ) {
          int j = D0.O(i);
          int j_inv = inv[j];
          int j_inv_pr = D0.O( dinv[i] );
          if ( j_inv_pr != j_inv ) {
               FatalErr( "dinv for D0 not consistent with inv for hb" );
          }
     }


#if 0
     // are there dinv[i] == i ???
     for ( int i = 0; i < dinv.isize(); ++i ) {
          if ( dinv[i] == i ) {
               cout << "dinv[" << i << "]=" << dinv[i] << endl;
               int j = D0.O(i);
               cout << "    j=" << j << ", inv[j]=" << inv[j] << endl;

               auto fw = hb.O(j);
               auto rc = fw;
               rc.ReverseComplement();
               if ( fw != rc ) {
                    cout << "fw != rc" << endl;
                    cout << "fw: ";
                    fw.PrintBases( cout, 0, 20);
                    cout << endl;
                    cout << "rc: ";
                    rc.PrintBases(cout, 0, 20);
               } else {
                    cout << "fw == rc OK!!!!" << endl;
               }
          }
     }
#endif

     // validate D0 connectivity
     vec<int> to_left, to_right;
     hb.ToLeft(to_left);
     hb.ToRight(to_right);

     cout << Date() << ": checking D0 consistency with hb" << endl;
     for ( int v = 0; v < D0.N(); ++v ) {
          for ( int i = 0; i < D0.To(v).isize(); ++i )
               for ( int j = 0; j < D0.From(v).isize(); ++j ) {
                    int to_edge = D0.OFrom(v, j);
                    int from_edge = D0.OTo(v, i);

                    if ( to_edge != from_edge &&
                              to_right[from_edge] != to_left[to_edge] ) {
                         PRINT3( v, to_edge, from_edge );
                         FatalErr( "D0 is inconsistent with hb");
                    }
               }

     }
}

void ValidateClosures( vec<vec<int>> const& all_closures, HyperBasevector const& hb, vec<int> const& inv )
{
     size_t bads = 0;

     // make an index of edge starts
     vec<vec<int>> ei( hb.E() );

     for ( int i = 0; i < all_closures.isize(); ++i ) {
          ei[all_closures[i][0]].push_back( i );
     }

     for ( auto& vi : ei ) { UniqueSort(vi); }

     // for each closure, find its rc
     for ( int i = 0; i < all_closures.isize(); ++i ) {
          auto const& fw = all_closures[i];
          vec<int> rc;
          for ( int j = fw.size() - 1; j >=0; j-- )
               rc.push_back( inv[fw[j]] );

          // now find the rc
          bool found = false;
          for ( size_t j = 0; j < ei[rc[0]].size() ; j++ ) {
               if ( rc == all_closures[ei[rc[0]][j]] ) {
                    found = true;
                    break;
               }
          }
          if ( !found ) bads++;
     }

     // for each closure, validate its path in hb
     vec<int> to_left, to_right;
     hb.ToLeft(to_left);
     hb.ToRight(to_right);

     size_t bads2 = 0;
     for ( int i = 0; i < all_closures.isize(); ++i ) {
          if ( !ValidatePath( all_closures[i], hb, to_left, to_right ) ) bads2++;
     }

     cout << Date() << ": Closure validation: " << all_closures.size() << " closures, "
          << bads << " missing rc, " << bads2 << " bad paths." << endl;
}

void MergeClosure( vec<vec<int>>& AC,
     vec<vec<om>>& omatch, const int c1, int c2, int start1, int start2, int len )
{
     // Merge.

     vec<int> &C1 = AC[c1], &C2 = AC[c2];
     vec<int> C1new;
     int left_ext1 = 0, left_ext2 = 0;
     if ( start1 >= start2 )
     {    for ( int j = 0; j < start1; j++ )
               C1new.push_back( C1[j] );
          left_ext2 = start1 - start2;    }
     else
     {    for ( int j = 0; j < start2; j++ )
               C1new.push_back( C2[j] );
          left_ext1 = start2 - start1;    }
     if ( C1.isize( ) - start1 >= C2.isize( ) - start2 )
     {    for ( int j = start1; j < C1.isize( ); j++ )
               C1new.push_back( C1[j] );    }
     else
     {    for ( int j = start2; j < C2.isize( ); j++ )
               C1new.push_back( C2[j] );    }
     C1 = C1new;
     if ( left_ext1 > 0 )
     {    for ( int j = 0; j < omatch[c1].isize( ); j++ )
               omatch[c1][j].start1 += left_ext1;
          vec<int> c3s;
          for ( int j = 0; j < omatch[c1].isize( ); j++ )
          {    int c3 = omatch[c1][j].c2;
               if ( AC[c3].nonempty( ) ) c3s.push_back(c3);    }
          UniqueSort(c3s);
          for ( int i = 0; i < c3s.isize( ); i++ )
          {    int c3 = c3s[i];
               for ( int j = 0; j < omatch[c3].isize( ); j++ )
               {    if ( omatch[c3][j].c2 != c1 ) continue;
                    omatch[c3][j].start2 += left_ext1;    }    }    }
     for ( int l = 0; l < omatch[c2].isize( ); l++ )
     {
          // Translate match.

          om& o2 = omatch[c2][l];
          int c3 = o2.c2;
          const vec<int>& C3 = AC[c3];
          if ( C3.empty( ) ) continue;
          int s1 = o2.start1 + left_ext2, s3 = o2.start2, len = o2.len;
          om x( c3, s1, s3, len );
          x.Extend( C1, C3 );

          // See if we already have it.

          Bool found = False;
          for ( int m = 0; m < omatch[c1].isize( ); m++ )
          {    if ( omatch[c1][m].c2 != c3 ) continue;
               omatch[c1][m].Extend( C1, C3 );
               if ( x == omatch[c1][m] )
               {    found = True;
                    break;    }    }
          if (found) continue;

          // Save new match.

          omatch[c1].push_back(x);
          omatch[c3].push( c1, s3, s1, len );    }

     // Clean up.

     C2.clear( );
     omatch[c2].clear( );    }

void ClosuresToGraph( const HyperBasevectorX& hb, const vec<int>& inv,
     vec<vec<int>>& all_closures, digraphE<vec<int>>& D, vec<int>& dinv,
     const Bool verbose )
{
     // Index closures.

     int64_t N = all_closures.size( ), total = 0;
     for ( int64_t i = 0; i < N; i++ )
          total += all_closures[i].size( );
     if (verbose)
     {    cout << Date( ) << ": " << ToStringAddCommas(N) << " closures "
               << "having in total " << ToStringAddCommas(total) << " edges" << endl;
          cout << Date( ) << ": indexing closures, peak mem = "
               << PeakMemUsageGBString( ) << endl;    }
     vec<vec<int>> ci( hb.E( ) );
     for ( int64_t i = 0; i < N; i++ )
     for ( int j = 0; j < all_closures[i].isize( ); j++ )
          ci[ all_closures[i][j] ].push_back(i);
     if (verbose) cout << Date( ) << ": uniquesorting" << endl;
     #pragma omp parallel for schedule( dynamic, 1000 )
     for ( int e = 0; e < hb.E( ); e++ )
          UniqueSort( ci[e] );

     // Set up match data structure.

     vec<vec<om>> omatch(N); // { (c2,start1,start2,len) }

     // Find matches.  These are defined by extensions of one closure by another.
     // There are other matches which could in principle be mined.

     if (verbose)
     {    cout << Date( ) << ": finding overlaps, peak mem = "
               << PeakMemUsageGBString( ) << endl;    }
     const int MIN_OVER = 200 - (hb.K()-1);
     #pragma omp parallel for schedule( dynamic, 1000 )
     for ( int i1 = 0; i1 < N; i1++ )
     {    const vec<int>& x1 = all_closures[i1];
          int nkmers = 0, b = -1, best = 1000000000;
          for ( int j = x1.isize( ) - 1; j >= 0; j-- )
          {    if ( ci[ x1[j] ].isize( ) < best )
               {    best = ci[ x1[j] ].size( );
                    b = j;    }
               nkmers += hb.Kmers( x1[j] );
               if ( nkmers >= MIN_OVER ) break;    }
          const vec<int>& ids = ci[ x1[b] ];
          for ( int u = 0; u < ids.isize( ); u++ )
          {    int i2 = ids[u];
               const vec<int>& x2 = all_closures[i2];
               int j1 = b;
               int stop2 = Min( x2.isize( ), j1 + x2.isize( ) - x1.isize( ) );
               for ( int j2 = 0; j2 < stop2; j2++ )
               {    if ( x1[j1] != x2[j2] ) continue;
                    Bool match = True;
                    int over = 0;
                    int start1 = Max( 0, j1 - j2 );
                    int stop1 = Min( x1.isize( ), x2.isize( ) + j1 - j2 );
                    for ( int i1 = start1; i1 < stop1; i1++ )
                    {    int i2 = j2 - (j1-i1);
                         if ( x1[i1] != x2[i2] )
                         {    match = False;
                              break;    }
                         over += hb.Kmers( x1[i1] );    }
                    if ( !match || over < MIN_OVER ) continue;
                    int len = stop1 - start1;
                    int start2 = j2 - (j1-start1);
                    omatch[i1].push( i2, start1, start2, len );    }    }    }

     // Add all long matches.

     if (verbose) cout << Date( ) << ": adding long matches" << endl;
     const int batch = 10000;
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int bi = 0; bi < hb.E( ); bi += batch )
     {    vec< pair<int,int> > Q;
          for ( int e = bi; e < Min( bi + batch, hb.E( ) ); e++ )
          {    if ( hb.Kmers(e) < MIN_OVER ) continue;
               if ( ci[e].size( ) <= 1 ) continue;

               // Locate e.

               Q.clear( );
               for ( int j = 0; j < ci[e].isize( ); j++ )
               {    int c = ci[e][j];
                    const vec<int>& x = all_closures[c];
                    for ( int m = 0; m < x.isize( ); m++ )
                         if ( x[m] == e ) Q.push( c, m );    }

               // Join.

               for ( int i1 = 0; i1 < Q.isize( ); i1++ )
               {    int c1 = Q[i1].first, m1 = Q[i1].second;
                    for ( int i2 = i1 + 1; i2 < Q.isize( ); i2++ )
                    {    int c2 = Q[i2].first, m2 = Q[i2].second;
                         Bool match = False;
                         for ( int l = 0; l < omatch[c1].isize( ); l++ )
                         {    if ( omatch[c1][l].c2 != c2 ) continue;
                              if ( omatch[c1][l].Offset( ) != m1 - m2 ) continue;
                              if ( m1 < omatch[c1][l].Start1( ) ) continue;
                              if ( m1 >= omatch[c1][l].Stop1( ) ) continue;
                              match = True;
                              break;    }
                         if (match) continue;
                         #pragma omp critical
                         {    omatch[c1].push( c2, m1, m2, 1 );
                              omatch[c1].back( ).Extend(
                                   all_closures[c1], all_closures[c2] );
                                   }    }    }    }    }

     // Force match symmetry.

     if (verbose) cout << Date( ) << ": symmetrizing" << endl;
     vec<int> oms(N);
     for ( int64_t i = 0; i < N; i++ )
          oms[i] = omatch[i].size( );
     for ( int64_t i1 = 0; i1 < N; i1++ )
     {    for ( int j = 0; j < oms[i1]; j++ )
          {    int64_t i2 = omatch[i1][j].c2;
               int start1 = omatch[i1][j].start1, start2 = omatch[i1][j].start2;
               int len = omatch[i1][j].len;
               omatch[i2].push( i1, start2, start1, len );    }    }

     // Compute involution of closures.

     if (verbose)
          cout << Date( ) << ": computing involution of all_closures" << endl;
     vec<int64_t> cinv( N, -1 );
     #pragma omp parallel for
     for ( int64_t i = 0; i < N; i++ )
     {    const vec<int>& x = all_closures[i];
          vec<int> rx( x.size( ) );
          for ( int j = 0; j < x.isize( ); j++ )
               rx[ x.isize( ) - j - 1 ] = inv[ x[j] ];
          int64_t ip = BinPosition( all_closures, rx );
          ForceAssertGe( ip, 0 );
          cinv[i] = ip;    }

     // Force match symmetry with respect to involution.

     for ( int64_t i = 0; i < N; i++ )
          oms[i] = omatch[i].size( );
     if (verbose) cout << Date( ) << ": finding adds" << endl;
     int64_t iadds = 0, itotal = 0;
     double iclock = WallClockTime( );
     for ( int64_t i1 = 0; i1 < N; i1++ )
     {    itotal += oms[i1];
          for ( int j = 0; j < oms[i1]; j++ )
          {    int64_t i2 = omatch[i1][j].c2;
               int start1 = omatch[i1][j].start1, start2 = omatch[i1][j].start2;
               int len = omatch[i1][j].len;
               int64_t ip1 = cinv[i1], ip2 = cinv[i2];
               int istart1 = all_closures[i1].isize( ) - start1 - len;
               int istart2 = all_closures[i2].isize( ) - start2 - len;
               om o( ip2, istart1, istart2, len );
               if ( !Member( omatch[ip1], o ) )
               {    omatch[ip1].push_back(o);
                    iadds++;    }    }    }
     if (verbose)
     {    DPRINT2_TO( cout, iadds, itotal );
          cout << Date( ) << ": done, time used = " << TimeSince(iclock)
               << endl;    }

     /*
     if (WRITE)
     {    if (verbose) cout << Date( ) << ": writing matches" << endl;
          BinaryWriter::writeFile( DIR + "/a.omatch" + WRITE_SUFFIX, omatch );    }
     */

     // Partially merge closures by modifying all_closures and omatch.  The purpose
     // of this step is solely to reduce memory usage by the full-blown gluing
     // step of the algorithm.

     vec<vec<int>>& AC = all_closures;
     if (verbose)
     {    cout << Date( ) << ": start merging closures, mem = "
               << MemUsageGBString( ) << endl;    }
     {    int64_t merges = 0;
          for ( int64_t c1 = 0; c1 < AC.jsize( ); c1++ )
          {    if ( AC[c1].empty( ) ) continue;
               vec<Bool> od1( omatch[c1].size( ), False );
               vec<om> add1;
               for ( int i = 0; i < omatch[c1].isize( ); i++ )
               {    int c2 = omatch[c1][i].c2;
                    if ( c2 == c1 ) continue;
                    if ( AC[c2].empty( ) ) continue;
                    int& start1 = omatch[c1][i].start1;
                    int& start2 = omatch[c1][i].start2;
                    int& len = omatch[c1][i].len;
                    const vec<int> &C1 = AC[c1], &C2 = AC[c2];
                    int n1 = C1.size( ), n2 = C2.size( );

                    // Extend the match, then test for proper overlap.

                    omatch[c1][i].Extend( C1, C2 );
                    if ( start1 > 0 && start2 > 0 ) continue;
                    if ( start1 + len < n1 && start2 + len < n2 ) continue;

                    // Require that c1 is to the left of c2.
                    //
                    //         start1
                    // ----------------c1--------------
                    //         start2
                    //         ----------------c2----------------

                    if ( start2 > 0 ) continue;

                    // Merge.  Force symmetry with respect to involution.

                    int64_t ic1 = cinv[c1], ic2 = cinv[c2];
                    vec<int> &IC1 = all_closures[ic1], &IC2 = all_closures[ic2];
                    int il = -1;
                    om& y = omatch[c1][i];
                    for ( int l = 0; l < omatch[ic1].isize( ); l++ )
                    {    const om& x = omatch[ic1][l];
                         if ( x.c2 != ic2 ) continue;
                         if ( x.Offset( ) !=
                              ( C1.isize( ) - y.Offset( ) ) - C2.isize( ) )
                         {    continue;    }
                         int istart1 = C1.isize( ) - y.start1 - y.len;
                         int istart2 = C2.isize( ) - y.start2 - y.len;
                         if ( IntervalOverlap( x.start1, x.start1 + x.len,
                              istart1, istart1 + y.len ) <= 0 )
                         {    continue;    }
                         il = l;
                         break;    }
                    ForceAssert( il >= 0 );
                    omatch[ic1][il].Extend( IC1, IC2 );
                    if ( cinv[c1] == c1 || cinv[c1] == c2
                         || cinv[c2] == c1 || cinv[c2] == c2 )
                    {    continue;    }
                    MergeClosure( AC, omatch, c1, c2, omatch[c1][i].start1,
                         omatch[c1][i].start2, omatch[c1][i].len );
                    MergeClosure( AC, omatch, ic1, ic2, omatch[ic1][il].start1,
                         omatch[ic1][il].start2, omatch[ic1][il].len );
                    merges += 2;    }    }

          // Clean up.

          if (verbose)
          {    cout << Date( ) << ": cleaning up merger, mem = "
                    << MemUsageGBString( ) << endl;    }
          for ( int64_t c1 = 0; c1 < AC.jsize( ); c1++ )
          {    if ( AC[c1].empty( ) ) omatch[c1].clear( );
               else
               {    vec<Bool> to_delete( omatch[c1].size( ), False );
                    for ( int j = 0; j < omatch[c1].isize( ); j++ )
                    {    int c2 = omatch[c1][j].c2;
                         if ( AC[c2].empty( ) ) to_delete[j] = True;    }
                    EraseIf( omatch[c1], to_delete );    }    }
          if (verbose)
          {    cout << Date( ) << ": merge made "
                    << ToStringAddCommas(merges) << " merges in total"
                    << endl;    }    }

     // Clean up.

     if (verbose)
     {    cout << Date( ) << ": start clean up after merger, mem = "
               << MemUsageGBString( ) << endl;    }
     {    vec<int64_t> to_new( AC.size( ), -1 );
          vec<Bool> cto_delete( AC.size( ), False );
          int64_t acc = 0;
          for ( int64_t i = 0; i < AC.isize( ); i++ )
          {    if ( AC[i].empty( ) ) cto_delete[i] = True;
               else
               {    to_new[i] = acc;
                    acc++;    }    }
          EraseIf( AC, cto_delete ), EraseIf( omatch, cto_delete );
          EraseIf( cinv, cto_delete );
          #pragma omp parallel for
          for ( int64_t i = 0; i < omatch.jsize( ); i++ )
          {    vec<Bool> to_delete( omatch[i].size( ), False );
               for ( int j = 0; j < omatch[i].isize( ); j++ )
               {    int& c2 = omatch[i][j].c2;
                    if ( to_new[c2] < 0 ) to_delete[j] = True;
                    else c2 = to_new[c2];    }    }
          for ( int64_t i = 0; i < cinv.jsize( ); i++ )
               cinv[i] = to_new[ cinv[i] ];    }

     // Extend

     if (verbose)
     {    cout << Date( ) << ": extending overlaps, mem = "
               << MemUsageGBString( ) << endl;    }
     for ( int64_t c1 = 0; c1 < AC.isize( ); c1++ )
     for ( int j = 0; j < omatch[c1].isize( ); j++ )
     {    int c2 = omatch[c1][j].c2;
          const vec<int>& C1 = AC[c1];
          const vec<int>& C2 = AC[c2];
          omatch[c1][j].Extend( C1, C2 );    }

     // Report stats.

     if (verbose)
     {    N = all_closures.size( ), total = 0;
          for ( int64_t i = 0; i < N; i++ )
               total += all_closures[i].size( );
          cout << Date( ) << ": " << ToStringAddCommas(N) << " closures "
               << "having in total " << ToStringAddCommas(total) << " edges" << endl;
          cout << Date( ) << ": done merging closures" << endl;    }

     // Merge words.

     if (verbose)
     {    cout << Date( ) << ": merging, peak mem = " << PeakMemUsageGBString( )
               << endl;    }
     digraphE<int>* D0p = new digraphE<int>;
     // ValidateClosures( all_closures, hb, inv );
     WDG( all_closures, cinv, omatch, *D0p, dinv, verbose );
     // ValidateD0( hb, inv, *D0p, dinv );
     Bool single = False;
     Vectorify( D0p, dinv, D, verbose, single, True );    }
