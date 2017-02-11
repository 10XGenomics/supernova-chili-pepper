// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Limitations.  This will not seed on sequence gap edges in D.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "PackAlign.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/HyperBasevector.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/large/Lines.h"
#include "random/Shuffle.h"
#include "10X/DfTools.h"
#include "10X/Heuristics.h"
#include "10X/astats/N50P.h"

void Align( const basevector& X1, const basevector& X2, int& score,
     vec<ho_interval>& perfs2 )
{    
     // Heuristics.

     const int MAX_DELTA = 40;
     const int BW_ADD = 5;
     const int MIN_TO_SUPER = 15000;

     // Align.

     alignment al;
     int best_loc, offset = X1.isize( ) - X2.isize( );
     align a;
     perfs2.clear( );
     if ( X1.size( ) == 0 || X2.size( ) == 0 ) 
     {    score = 3 * Max( X1.size( ), X2.size( ) );    }
     else if ( Abs(offset) <= MAX_DELTA )
     {    int bandwidth = Abs(offset) + BW_ADD;
          int errors;
          score = 2 * SmithWatBandedA( X1, X2, offset, bandwidth, a, errors );
          a.PerfectIntervals2( X1, X2, perfs2 );    }
     else if ( X2.size( ) >= MIN_TO_SUPER )
     {    const int64_t max_product = 1000000000;
          uint answer = SmithWatAffineSuper( X1, X2, al, 
               // next line = copying default args
               501, true, true, 0, 3, 12, 1, SmithWatAffineParallel2, 
               max_product );
          if ( answer != SmithWatAffineSuper_FAIL ) 
          {    score = answer;
               a = al;
               a.PerfectIntervals2( X1, X2, perfs2 );    }
          else score = 1000000000;    }
     else if ( X2.size( ) > X1.size( ) )
     {    score = SmithWatFree( X1, X2, best_loc, al, true, true );
          a = al;
          a.PerfectIntervals2( X1, X2, perfs2 );    }
     else
     {    score = SmithWatFree( X2, X1, best_loc, al, true, true );    
          a = al;
          a.PerfectIntervals1( X2, X1, perfs2 );    }    }

int64_t N50P( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<int>& dinv, const vecbasevector& G, 
     const vec<vec<pair<int,int>>>& locs, const vec<Bool>& keep, 
     String& report, String& details, double& errw )
{
     // Translate locs to this graph.
     //
     // hlocs[g][ index in locs[g] ]  = ( edge in D, start on D )

     vec<vec<vec<pair<int,int>>>> hlocs( locs.size( ) );
     cout << Date( ) << ": pushing" << endl;
     {    vec<vec<pair<int,int>>> hpos( hb.E( ) );
          for ( int d = 0; d < D.E( ); d++ )
          {    if ( !keep[d] ) continue;
               const vec<int>& x = D.O(d);
               if ( x[0] < 0 ) continue;
               int n = 0;
               for ( int j = 0; j < x.isize( ); j++ )
               {    hpos[ x[j] ].push( d, n );
                    n += hb.Kmers( x[j] );    }    }
          cout << Date( ) << ": translating locs" << endl;
          #pragma omp parallel for
          for ( int g = 0; g < locs.isize( ); g++ )
          {    hlocs[g].resize( locs[g].size( ) );
               for ( int p = 0; p < locs[g].isize( ); p++ )
               {    int e = locs[g][p].first, epos = locs[g][p].second;
                    if ( e < 0 ) continue;
                    for ( int j = 0; j < hpos[e].isize( ); j++ )
                    {    int pos = epos + hpos[e][j].second;
                         int d = hpos[e][j].first;
                         hlocs[g][p].push( d, pos );    } 
                    UniqueSort( hlocs[g][p] );    }    }    }

     // Find lines.

     cout << Date( ) << ": finding lines" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     vec<pair<int,int>> tol( D.E( ), make_pair(-1,-1) );
     for ( int i = 0; i < dlines.isize( ); i++ )
     for ( int j = 0; j < dlines[i].isize( ); j++ )
     for ( int k = 0; k < dlines[i][j].isize( ); k++ )
     for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
          tol[ dlines[i][j][k][l] ] = make_pair( i, j );

     // Convert to HyperBasevector.

     cout << Date( ) << ": converting to hbd" << endl;
     HyperBasevectorX hbd;
     SuperToSeqGraph( hb, D, hbd );

     // Create a set of match blocks.  This is not done completely correctly.
     // Perhaps better, the match blocks should be required to satisfy a 
     // unipath-like condition, going up to but not over branches.  Note that 
     // branching can only arise from failures of zippering.  These could arise at 
     // sequence "gaps" and in rare instances otherwise.

     cout << Date( ) << ": finding match blocks" << endl;
     const int IGNORE = 300;
     vec<vec< quad<ho_interval,vec<int>,int,int> >> PERFS( locs.size( ) );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    vec<vec<Bool>> used( hlocs[g].size( ) );
          for ( int i = 0; i < hlocs[g].isize( ); i++ )
               used[i].resize( hlocs[g][i].size( ), False );
          for ( int start = 0; start < locs[g].isize( ); start++ )
          for ( int j = 0; j < hlocs[g][start].isize( ); j++ )
          {    if ( used[start][j] ) continue;
               int stop;
               int h = hlocs[g][start][j].first;       // edge in D
               int pos = hlocs[g][start][j].second;    // start on D
               int start_pos = pos;
               vec<int> p = {h};
               // test here and below needed because of the trimming
               // induced by seq gap edges
               if ( pos + hb.K( ) > hbd.Bases(h) ) continue;
               for ( stop = start+1; stop < locs[g].isize( ); stop++ )
               {    int K = hb.K( );

                    if ( pos + K < hbd.Bases(h) )
                    {    if ( G[g][stop+K-1] == hbd.O(h)[pos+K] )
                         {    int k = BinPosition( 
                                   hlocs[g][stop], make_pair( h, pos + 1 ) );
                              if ( k >= 0 )
                              {    pos++;
                                   used[stop][k] = True;
                                   continue;    }    }    }

                    Bool ext = False;
                    for ( int k = 0; k < hlocs[g][stop].isize( ); k++ )
                    {    if ( used[stop][k] ) continue;
                         int hp = hlocs[g][stop][k].first;
                         int posp = hlocs[g][stop][k].second;
                         if ( posp + hb.K( ) > hbd.Bases(hp) ) continue;
                         if ( hp == h && posp == pos + 1 )
                         {    pos++;
                              used[stop][k] = True;
                              ext = True;
                              break;    }
                         int v = hbd.ToRight(h), w = hbd.ToLeft(hp);
                         if ( v == w && pos == hbd.Bases(h) - hb.K( ) && posp == 0 )
                         {    h = hp;
                              p.push_back(hp);
                              pos = posp;
                              ext = True;
                              break;    }    }
                    if ( !ext ) break;    }
               int stop_pos = pos;
               int K = hb.K( );
               if ( stop - start > IGNORE )
               {    PERFS[g].push( ho_interval(start,stop+K-1), 
                         p, start_pos, stop_pos+K-1 );    }    }    }

     // Show matches.

     cout << Date( ) << ": printing match blocks" << endl;
     ostringstream dout;
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    dout << "\nMATCHES FOR " << g << endl;
          for ( int j = 0; j < PERFS[g].isize( ); j++ )
          {    dout << PERFS[g][j].first << ": " 
                    << printSeq( PERFS[g][j].second ) << " (" << PERFS[g][j].third 
                    << "," << PERFS[g][j].fourth << ")" << endl;    }    }

     // Merge match blocks.

     cout << Date( ) << ": merging match blocks for " << G.size( ) << " ref seqs" 
          << endl;
     vec<vec< quad<ho_interval,vec<int>,int,int> >> GLOGS = PERFS;
     int count = 0;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    const int MAX_PATHS = 10;
          const int MAX_ITERATIONS = 100;
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec< quad<ho_interval,vec<int>,int,int> > GLOGS2;
               int n = GLOGS[g].size( );
               vec<Bool> used( n, False );
               for ( int j = 0; j < n; j++ )
               {    if ( used[j] ) continue;
                    int start = GLOGS[g][j].first.Start( );
                    int stop = GLOGS[g][j].first.Stop( );
                    vec<int> x = GLOGS[g][j].second;
                    int left_stop = GLOGS[g][j].fourth;
                    int k;
                    for ( k = j + 1; k < n; k++ )
                    {    if ( used[k] ) continue;
                         vec<int> y = GLOGS[g][k].second;
                         if ( x.back( ) == y.front( ) )
                         {    stop = GLOGS[g][k].first.Stop( );
                              for ( int l = 1; l < y.isize( ); l++ )
                                   x.push_back( y[l] );
                              left_stop = GLOGS[g][k].fourth;
                              used[k] = True;    }
                         else if ( hbd.ToRight( x.back( ) ) 
                              == hbd.ToLeft( y.front( ) ) )
                         {    stop = GLOGS[g][k].first.Stop( );
                              x.append(y);
                              left_stop = GLOGS[g][k].fourth;
                              used[k] = True;    }
                         else
                         {    int right_start = GLOGS[g][k].third;
                              int xstart = stop, xstop = GLOGS[g][k].first.Start( );
     
                              // Data at this point:
                              // x = path we're using
                              // left_stop = position on last edge of x
                              // xstart = corresponding pos on ref 
                              // y = next path
                              // right_start = position on first edge of y
                              // xstop = corresponding pos on ref
     
                              int e1 = x.back( ), e2 = y.front( );
     
                              // If e1 and e2 are in the same cell, try to trim back
                              // to the boundaries.
     
                              int K = hbd.K( );
                              if ( tol[e1] == tol[e2] )
                              {    pair<int,int> p = tol[e1];
                                   while ( tol[e1] == p )
                                   {    if ( x.solo( ) ) break;
                                        int sub = left_stop - (K-2);
                                        Bool mismatch = False;
                                        for ( int l = 1; l < sub; l++ )
                                        {    if ( xstart - l < 0 || left_stop-l < 0
                                                  || G[g][xstart-l] 
                                                       != hbd.O(e1)[left_stop-l] )
                                             {    mismatch = True;
                                                  break;    }    }
                                        if (mismatch) break;
                                        x.pop_back( );
                                        xstart -= sub;
                                        e1 = x.back( );
                                        left_stop = hbd.Bases(e1) - 1;    }
                                   while ( tol[e2] == p )
                                   {    if ( y.size( ) == 1 ) break;
                                        Bool mismatch = False;

                                        int add 
                                             = hbd.Bases(e2) - (K-1) - right_start;
                                        if ( add < 0 )
                                        {    cout << "FUNNY!" << endl;
                                             break;    }
                                        for ( int l = 0; l < add; l++ )
                                        {    if ( xstop + l >= G[g].isize( )
                                                  || right_start + l >= hbd.Bases(e2)
                                                  || G[g][xstop+l]
                                                       != hbd.O(e2)[right_start+l] )
                                             {    mismatch = True;
                                                  break;    }    }
     
                                        // pos 0 on e2_next
                                        // = pos ( |e2| - (K-1) ) on e2
     
                                        if (mismatch) break;
                                        y.pop_front( );
                                        right_start = 0;
                                        xstop += add;
                                        e2 = y[0];    }    }

                              // Look for paths.
                    
                              xstop += K-1;
                              if ( xstop <= xstart ) continue;
                              int v = hbd.ToRight(e1), w = hbd.ToLeft(e2);
                              vec<vec<int>> paths;
                              Bool ok = hbd.EdgePaths( 
                                   v, w, paths, -1, MAX_PATHS, MAX_ITERATIONS );
                              if ( !ok || paths.empty( ) ) continue;
                              Bool gap = False;
                              for ( int i = 0; i < paths.isize( ); i++ )
                              for ( int l = 0; l < paths[i].isize( ); l++ )
                                   if ( hbd.Bases( paths[i][l] ) == 0 ) gap = True;
                              if (gap) continue;
     
                              vec<int> errs( paths.size( ), 0 );
                              right_start += K-1;

                              basevector GG( G[g], xstart, xstop - xstart );
                              vecbasevector bpaths;
                              for ( int l = 0; l < paths.isize( ); l++ )
                              {    basevector b = hbd.Cat( paths[l] );
                                   b.SetToSubOf( b, K-1, b.isize( ) - (K-1) );
                                   bpaths.push_back(b);    }
                              if ( left_stop < hbd.Bases(e1) )
                              {    basevector left( hbd.O(e1), left_stop, 
                                        hbd.Bases(e1) - left_stop );
                                   for ( int l = 0; l < paths.isize( ); l++ )
                                   {    basevector b = left;
                                        b.append( bpaths[l] );
                                        bpaths[l] = b;    }    }
                              if ( right_start > K-1 )
                              {    basevector right(
                                              hbd.O(e2), K-1, right_start - (K-1) );
                                   for ( int l = 0; l < paths.isize( ); l++ )
                                        bpaths[l].append(right);    }
                              else if ( right_start < K-1 )
                              {    for ( int l = 0; l < paths.isize( ); l++ )
                                   {    bpaths[l].resize( bpaths[l].isize( ) 
                                             - ( (K-1) - right_start ) );    }    }
                              vec<int> score( paths.size( ) );
                              for ( int l = 0; l < paths.isize( ); l++ )
                              {    vec<ho_interval> perfs2;
                                   Align( bpaths[l], GG, score[l], perfs2 );    }
                              SortSync( score, paths );
                              x.append( paths[0] );
                              for ( auto e : y ) x.push_back(e);
                              stop = GLOGS[g][k].first.Stop( );
                              left_stop = GLOGS[g][k].fourth;    
                              used[k] = True;    }    }
                    GLOGS2.push( ho_interval(start,stop), x, 
                         GLOGS[g][j].third, left_stop );    }
               GLOGS[g] = GLOGS2;    }

          // Update status.

          #pragma omp critical
          {    if ( count > 0 && count % 50 == 0 ) cout << "\n";
               else if ( count > 0 && count % 10 == 0 ) cout << " ";
               cout << ".";
               flush(cout);
               count++;    }    }
     cout << endl;

     // Show condensed matches.

     cout << Date( ) << ": printing condensed match blocks" << endl;
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    dout << "\nCONDENSED MATCHES FOR " << g << endl;
          for ( int j = 0; j < GLOGS[g].isize( ); j++ )
          {    dout << GLOGS[g][j].first << ": " 
                    << printSeq( GLOGS[g][j].second ) << " (" << GLOGS[g][j].third 
                    << "," << GLOGS[g][j].fourth << ")" << endl;    }    }

     // Now combine across gaps.

     cout << Date( ) << ": combining across gaps" << endl;
     vec<vec<vec< quad<ho_interval,vec<int>,int,int> >>> soos( G.size( ) );
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    vec<Bool> used( GLOGS[g].size( ), False );
          for ( int j = 0; j < GLOGS[g].isize( ); j++ )
          {    if ( used[j] ) continue;
               vec< quad<ho_interval,vec<int>,int,int> > S = { GLOGS[g][j] };
               for ( int k = j + 1; k < GLOGS[g].isize( ); k++ )
               {    const vec<int>& p1 = S.back( ).second, p2 = GLOGS[g][k].second;
                    int d1 = p1.back( ), d2 = p2.front( );
                    int v = hbd.ToRight(d1), w = hbd.ToLeft(d2);
                    if ( hbd.From(v).size( ) != 1 ) continue;
                    if ( hbd.To(w).size( ) != 1 ) continue;
                    if ( hbd.From(v)[0] != w ) continue;
                    int d = hbd.IFrom(v,0);
                    if ( hbd.Bases(d) != 0 ) continue;
                    used[k] = True;
                    S.push_back( GLOGS[g][k] );    }
               soos[g].push_back(S);    }    }

     // Allow more complex gaps.
     // ---left---> ---gap---> ---edge---> ---gap---> ---right--->

     vec<vec<vec< quad<ho_interval,vec<int>,int,int> >>> soos2( G.size( ) );
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    vec<Bool> used( soos[g].size( ), False );
          for ( int j = 0; j < soos[g].isize( ); j++ )
          {    if ( used[j] ) continue;
               vec< quad<ho_interval,vec<int>,int,int> > S = soos[g][j];
               for ( int k = j + 1; k < soos[g].isize( ); k++ )
               {    const vec<int>& p1 = S.back( ).second; 
                    const vec<int>& p2 = soos[g][k].front( ).second;
                    int d1 = p1.back( ), d2 = p2.front( );
                    int v = hbd.ToRight(d1), w = hbd.ToLeft(d2);
                    if ( hbd.From(v).size( ) != 1 ) continue;
                    if ( hbd.To(w).size( ) != 1 ) continue;
                    if ( hbd.OFrom(v,0).size( ) > 0 ) continue;
                    if ( hbd.OTo(w,0).size( ) > 0 ) continue;
                    int x = hbd.From(v)[0], y = hbd.To(w)[0];
                    if ( hbd.From(x).size( ) != 1 ) continue;
                    if ( hbd.To(y).size( ) != 1 ) continue;
                    if ( hbd.From(x)[0] != y ) continue;
                    used[k] = True;
                    S.append( soos[g][k] );    }
               soos2[g].push_back(S);    }    }
     soos = soos2;

     // Delete tiny combined blocks.

     for ( int g = 0; g < (int) G.size( ); g++ )
     {    vec<Bool> to_delete( soos[g].size( ), False );
          for ( int j = 0; j < soos[g].isize( ); j++ )
          {    int start = soos[g][j].front( ).first.Start( );
               int stop = soos[g][j].back( ).first.Stop( );
               if ( stop - start < 400 ) to_delete[j] = True;    }
          EraseIf( soos[g], to_delete );    }

     // Show combined condensed matches, and score them.

     cout << Date( ) << ": printing combined condensed match blocks, mem = " 
          << MemUsageGBString( ) << endl;
     vec<vec<int>> scores( G.size( ) );
     vec<vec<vec<ho_interval>>> P( G.size( ) );
     vec<vec<Bool>> align_fail( G.size( ) );
     vec<String> creports( G.size( ) );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    ostringstream mout;
          mout << "\nCOMBINED CONDENSED MATCHES FOR " << g 
               << " (L=" << G[g].size( ) << ")" << endl;
          scores[g].resize( soos[g].size( ) );
          P[g].resize( soos[g].size( ) );
          align_fail[g].resize( soos[g].size( ), False );
          for ( int j = 0; j < soos[g].isize( ); j++ )
          {    mout << "[" << j+1 << "] ";
               mout << soos[g][j].front( ).first.Start( ) << "-"
                    << soos[g][j].back( ).first.Stop( ) << " = ";
               for ( int i = 0; i < soos[g][j].isize( ); i++ )
               {    if ( i > 0 ) mout << "; ";
                    mout << soos[g][j][i].first << ": " 
                         << printSeq( soos[g][j][i].second ) 
                         << " (" << soos[g][j][i].third << "," 
                         << soos[g][j][i].fourth << ")";    } 
               mout << endl;

               // Score the match.

               vec<ho_interval> perfs;
               const vec< quad<ho_interval,vec<int>,int,int> >& glogs = soos[g][j];
               for ( int u = 0; u < glogs.isize( ); u++ )
               {    int start = glogs[u].first.Start( ); 
                    int stop = glogs[u].first.Stop( );
                    vec<int> x = glogs[u].second;
                    basevector b = hbd.Cat(x);
                    // Not sure why we're adding one to astop in this line:
                    int astart = glogs[u].third, astop = glogs[u].fourth + 1;
                    int rtrim = hbd.Bases( x.back( ) ) - astop;
                    // Not sure how this could happen:
                    if ( rtrim < 0 || b.isize( ) - astart - rtrim < 0 )
                    {    align_fail[g][j] = True;
                         break;    }
                    b.SetToSubOf( b, astart, b.isize( ) - astart - rtrim );
                    basevector GG( G[g], start, stop - start );

                    // Try to avoid actual alignment.

                    Bool aligned = False;
                    if ( b.size( ) == GG.size( ) )
                    {    int min_dist = 1000000000, last_diff = -1;
                         vec<ho_interval> P;
                         for ( int j = 0; j < b.isize( ); j++ )
                         {    if ( GG[j] != b[j] )
                              {    if ( last_diff < 0 ) 
                                   {    if ( j > 0 ) P.push( start, start + j );    }
                                   else if ( last_diff + 1 < j )
                                   {    P.push( 
                                             start + last_diff + 1, start + j );    }
                                   min_dist = Min( min_dist, j - last_diff );
                                   last_diff = j;    }    }
                         if ( last_diff + 1 < b.isize( ) )
                              P.push( start + last_diff + 1, start + b.isize( ) );
                         if ( min_dist >= 20 )
                         {    perfs.append(P);
                              aligned = True;    }    }

                    // Align.

                    if ( !aligned )
                    {    vec<ho_interval> perfsu;
                         int score;
                         Align( b, GG, score, perfsu );
                         if ( score == 1000000000 )
                         {    align_fail[g][j] = True;
                              break;    }
                         for ( int i = 0; i < perfsu.isize( ); i++ )
                              perfsu[i].Shift(start);
                         perfs.append(perfsu);    }    }
               if ( align_fail[g][j] )
               {    mout << "Alignment failed." << endl;
                    continue;    }

               // Remove overlaps between perfect match intervals.

               Sort(perfs);
               vec<Bool> to_delete( perfs.size( ), False );
               for ( int i = 0; i < perfs.isize( ) - 1; i++ )
               {    if ( perfs[i].Stop( ) > perfs[i+1].Start( ) )
                    {    perfs[i].AddToStop( 
                              perfs[i+1].Start( ) - perfs[i].Stop( ) );    
                         if ( perfs[i].Length( ) == 0 ) 
                              to_delete[i] = True;    }    }
               EraseIf( perfs, to_delete );

               // Add missing points.

               if ( perfs.empty( ) )
               {    for ( int j = 0; j < G[g].isize( ); j++ )
                         perfs.push( j, j+1 );    }
               else
               {    int np = perfs.size( );
                    for ( int j = 0; j < perfs[0].Start( ); j++ )
                         perfs.push( j, j+1 );
                    for ( int i = 0; i < np - 1; i++ )
                    for ( int j = perfs[i].Stop( ); j < perfs[i+1].Start( ); j++ )
                         perfs.push( j, j+1 );
                    for ( int j = perfs[np-1].Stop( ); j < G[g].isize( ); j++ )
                         perfs.push( j, j+1 );    }

               // Report.

               Sort(perfs);
               vec<int> perfsx;
               for ( int i = 0; i < perfs.isize( ); i++ )
                    perfsx.push_back( perfs[i].Length( ) );
               ReverseSort(perfsx);
               mout << "perfs = ";
               for ( int i = 0; i < perfsx.isize( ); i++ )
               {    if ( i > 0 ) mout << ",";
                    int j = perfsx.NextDiff(i);
                    if ( j - i == 1 ) mout << perfsx[i];
                    else mout << perfsx[i] << "[" << j-i << "]";
                    i = j - 1;    }
               scores[g][j] = N50(perfsx);
               P[g][j] = perfs;
               mout << "; N50 = " << scores[g][j] << endl;    }
          creports[g] = mout.str( );    }
     for ( int g = 0; g < (int) G.size( ); g++ )
          dout << creports[g];
     dout << endl;
     details = dout.str( );

     // Pick best match.

     cout << Date( ) << ": scoring matches, mem usage = " 
          << MemUsageGBString( ) << endl;
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    EraseIf( scores[g], align_fail[g] );
          EraseIf( soos[g], align_fail[g] );
          EraseIf( P[g], align_fail[g] );
          ReverseSortSync( scores[g], soos[g], P[g] );
          if ( soos[g].nonempty( ) ) 
          {    soos[g].resize(1);
               P[g].resize(1);    }    }

     // Handle the case where there's no alignment.

     for ( int g = 0; g < (int) G.size( ); g++ )
     {    if ( P[g].empty( ) )
          {    P[g].resize(1);
               for ( int j = 0; j < G[g].isize( ); j++ )
                    P[g][0].push( j, j+1 );    }    }

     // Build visual alignments.

     cout << Date( ) << ": building visuals, mem usage = " 
          << MemUsageGBString( ) << endl;
     vec<String> reports( G.size( ) );
     report += "\n";
     vec<pair<double,int>> err( G.size( ), make_pair(0,0) );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    ostringstream out;
          if ( g > 0 )
          {    out << "=========================================================="
                    << "==========================\n\n";    }
          out << "LOCATIONS FOR FINISHED SEQUENCE " << g 
               << " (LEN=" << G[g].size( ) << ")" << endl;
          if ( soos[g].empty( ) )
          {    out << "\n(none found)\n";
               continue;    }
          const vec< quad<ho_interval,vec<int>,int,int> >& glogs = soos[g][0];
          for ( int j = 0; j < glogs.isize( ); j++ )
          {    int start = glogs[j].first.Start( ), stop = glogs[j].first.Stop( );
               vec<int> x = glogs[j].second;
               out << endl << start << "-" << stop << ": " << printSeq(x) << endl;    
               basevector b = hbd.Cat(x);
               // Not sure why we're adding one to astop in this line:
               int astart = glogs[j].third, astop = glogs[j].fourth + 1;
               int rtrim = hbd.Bases( x.back( ) ) - astop;
               b.SetToSubOf( b, astart, b.isize( ) - astart - rtrim );
               basevector GG( G[g], start, stop - start );
               align a;
     
               // Try to avoid actual alignment.

               Bool aligned = False;
               if ( b.size( ) == GG.size( ) )
               {    int min_dist = 1000000000, last_diff = -1000000000;
                    for ( int j = 0; j < b.isize( ); j++ )
                    {    if ( GG[j] != b[j] )
                         {    min_dist = Min( min_dist, j - last_diff );
                              last_diff = j;    }    }
                    if ( min_dist >= 20 )
                    {    avector<int> gaps(1), lengths(1);
                         gaps(0) = 0;
                         lengths(0) = b.isize( );
                         a.Set( 0, 0, gaps, lengths );
                         aligned = True;    }    }

               // Align.
     
               if ( !aligned )
               {    alignment al;
                    SmithWatAffineSuper( b, GG, al );
                    a = al;    }

               // Print alignment.

               ostringstream pout;
               PrintVisualAlignmentClean( True, pout, b, GG, a );

               // Compute error penalty.

               const double mis_penalty = 1.0;
               const double ind_penalty = 1.0;
               const double ext_penalty = 0.1;
               int p1 = a.pos1( ), p2 = a.pos2( );
               for ( int j = 0; j < a.Nblocks( ); j++ ) 
               {    if ( a.Gaps(j) > 0 )  
                    {    err[g].first += ind_penalty;
                         err[g].first += ( a.Gaps(j) - 1 ) * ext_penalty;
                         err[g].second += a.Gaps(j);
                         p2 += a.Gaps(j);    }
                    if ( a.Gaps(j) < 0 )
                    {    err[g].first += ind_penalty;
                         err[g].first += ( -a.Gaps(j) - 1 ) * ext_penalty;
                         p1 -= a.Gaps(j);    }
                    err[g].second += a.Lengths(j);
                    for ( int x = 0; x < a.Lengths(j); x++ ) 
                    {    if ( b[p1] != GG[p2] ) err[g].first += mis_penalty;
                         ++p1;
                         ++p2;    }     }

               // Proceed.

               String s = pout.str( );
               if ( s.Contains( "(perfect", 0 ) ) s = "\n" + s;
               out << "\ntotal errs = " << err[g].first << "\n" << s;
               if ( j < glogs.isize( ) - 1 )
               {    int v = hbd.ToRight( x.back( ) );
                    if ( hbd.From(v).size( ) == 1 )
                    {    int e = hbd.IFrom(v,0);
                         if ( hbd.Bases(e) == 0 ) 
                         {    Bool two_sided = False;
                              const vec<int>& y = glogs[j+1].second;
                              int w = hbd.ToLeft( y.front( ) );
                              if ( hbd.To(w).size( ) == 1 )
                              {    int f = hbd.ITo(w,0);
                                   if ( f == e ) two_sided = True;    }
                              if (two_sided) out << "(two-sided gap)" << endl;
                              else out << "(gap, not sure to where)"
                                        << endl;    }    }    }    }
          reports[g] = out.str( );    }
     for ( int g = 0; g < (int) G.size( ); g++ )
          report += reports[g];

     // Compute weighted error.

     double num = 0, den = 0;
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    num += err[g].first;
          den += err[g].second;     }
     if ( den == 0 ) errw = 0;
     else errw = num/den;

     // Compute N50, allowing for random changing of the Fosmids.

     const int npasses = 1000;
     vec<int> N50s(npasses);
     vec<vec<int>> R(npasses);
     for ( int pass = 0; pass < npasses; pass++ )
     {    R[pass] = vec<int>( G.size( ), vec<int>::IDENTITY );
          Shuffle( R[pass].begin( ), R[pass].end( ) );    }
     cout << Date( ) << ": computing N50" << endl;
     #pragma omp parallel for
     for ( int pass = 0; pass < npasses; pass++ )
     {    vec<int> q;
          for ( int gg = 0; gg < (int) G.size( ); gg++ )
          {    int g = R[pass][gg];
               if ( gg > 0 )
               {    q.back( ) += P[g][0].front( ).Length( );
                    for ( int j = 1; j < P[g][0].isize( ); j++ )
                         q.push_back( P[g][0][j].Length( ) );    }
               else
               {    for ( int j = 0; j < P[g][0].isize( ); j++ )
                         q.push_back( P[g][0][j].Length( ) );    }    }
          N50s[pass] = N50(q);    }
     return int(round( Mean(N50s)));    }
