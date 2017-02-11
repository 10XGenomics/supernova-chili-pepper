// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "graph/DigraphTemplate.h"
#include "10X/WDG.h"

void MakeEquiv( const vec< vec<int> >& w, const vec<vec<om>>& omatch,
     const vec<int>& vstart, const vec<int>& estart,
     equiv_rel& evert, equiv_rel& eedge )
{
     for ( int i1 = 0; i1 < omatch.isize( ); i1++ )
     for ( int j = 0; j < omatch[i1].isize( ); j++ )
     {    int i2 = omatch[i1][j].c2;
          int start1 = omatch[i1][j].start1, start2 = omatch[i1][j].start2;
          int len = omatch[i1][j].len;
          for ( int l = 0; l <= len; l++ )
          {    int l1 = start1 + l, l2 = start2 + l;
               evert.Join( vstart[i1] + l1, vstart[i2] + l2 );
               if ( l < len )
                    eedge.Join( estart[i1] + l1, estart[i2] + l2 );    }    }    }

void WDG( vec< vec<int> >& w, vec<int64_t>& cinv, vec<vec<om>>& omatch, 
     digraphE<int>& D, vec<int>& dinv, const Bool verbose )
{    
     // Form the digraph W that is the disjoint union of the words, with each word
     // stretched out over a sequence of edges (corresponding to the letters in
     // the word).  Except that, we don't actually create W!  It exists only
     // virtually.

     double clock = WallClockTime( );
     vec< pair<int,int> > eorigin;
     vec<int> vstart, estart;
     int64_t iv, ie;
     if (verbose)
     {    cout << Date( ) << ": forming disjoint union, peak mem = " 
               << PeakMemUsageGBString( ) << endl;    }
     iv = 0, ie = 0;
     vstart.resize( w.size( ) );
     estart.resize( w.size( ) );
     int64_t nv = 0, ne = 0;
     for ( int i = 0; i < w.isize( ); i++ )
     {    ne += w[i].size( );
          nv += w[i].size( ) + 1;    }
     for ( int i = 0; i < w.isize( ); i++ )
     {    vstart[i] = iv, estart[i] = ie;
          for ( int j = 0; j < w[i].isize( ); j++ )
          {    eorigin.push( i, j );
               ie++;
               iv++;    }
          iv++;    }
     vec<int> ws( w.size( ) );
     #pragma omp parallel for
     for ( int64_t i = 0; i < w.jsize( ); i++ )
          ws[i] = w[i].size( );

     // From the overlaps, deduce equivalence relations on the vertices and edges
     // of W.

     equiv_rel evert(iv), eedge(ie);
     if (verbose)
     {    cout << Date( ) << ": deduce equivalence relation, peak mem = " 
               << PeakMemUsageGBString( ) << endl;    }
     MakeEquiv( w, omatch, vstart, estart, evert, eedge );
     Destroy(omatch);

     // Apply the equivalence relation to W, yielding the digraph D.

     if (verbose)
     {    cout << Date( ) << ": get orbit reps, peak mem = " 
               << PeakMemUsageGBString( ) << endl;    }
     vec<int> vreps, ereps;
     evert.OrbitRepsAlt(vreps), eedge.OrbitRepsAlt(ereps);
     if (verbose) 
     {    cout << Date( ) << ": form edits, peak mem = " 
               << PeakMemUsageGBString( ) << endl;    }
     int nvreps = vreps.size( ), nereps = ereps.size( );
     if (verbose)
     {    cout << Date( ) << ": setting edges2, peak mem = " 
               << PeakMemUsageGBString( ) << endl;    }
     D.EdgesMutable( ).resize(nereps);
     #pragma omp parallel for
     for ( int i = 0; i < nereps; i++ )
     {    int x = w[ eorigin[ ereps[i] ].first ][ eorigin[ ereps[i] ].second ];
          D.OMutable(i) = x;    }
     Destroy(w);
     D.FromMutable( ).resize(nvreps);
     D.ToMutable( ).resize(nvreps);
     D.FromEdgeObjMutable( ).resize(nvreps);
     D.ToEdgeObjMutable( ).resize(nvreps);

     // Create nex and en.  

     {    if (verbose) cout << Date( ) << ": creating en" << endl;
          int64_t nex = 0;
          vec<int> en(iv);
          vec<int64_t> startv( ws.size( ) ), starte( ws.size( ) );
          {    int64_t iv = 0, ie = 0;;
               for ( int i = 0; i < ws.isize( ); i++ )
               {    startv[i] = iv;
                    starte[i] = ie;
                    for ( int j = 0; j < ws[i]; j++ )
                    {    en[iv] = nex++;
                         ie++;
                         iv++;    }
                    iv++;    }    }

          // Apply equivalence relation.

          dinv.resize( ereps.size( ), -1 );
          if (verbose)
          {    cout << Date( ) << ": start || loop, peak mem = " 
                    << PeakMemUsageGBString( ) << endl;    }
          {    vec< triple<int,int,int> > edits(nex);
               vec< pair<int,int> > dinv_ed(nex);
               {    
                    #pragma omp parallel for schedule( dynamic, 10000 )
                    for ( int i = 0; i < ws.isize( ); i++ )
                    {    int64_t iv = startv[i], ie = starte[i];
                         for ( int j = 0; j < ws[i]; j++ )
                         {    int v1 = iv, v2 = iv+1, e = ie;
                              int pv1 = BinPosition( vreps, evert.ClassId(v1) );
                              int pv2 = BinPosition( vreps, evert.ClassId(v2) );
                              int pe = BinPosition( ereps, eedge.ClassId(e) );
                              edits[ en[v1] ] = make_triple( pv1, pv2, pe );

                              // Update graph inversion.

                              int64_t ri = cinv[i];
                              int rj = ws[i] - j - 1;
                              int re = estart[ri] + rj;
                              int rpe = BinPosition( ereps, eedge.ClassId(re) );
                              dinv_ed[ en[v1] ] = make_pair( pe, rpe );

                              // Advance.

                              ie++;
                              iv++;    }
                         iv++;    }    }

               // Apply inversion edits.

               for ( int64_t i = 0; i < nex; i++ )
                    dinv[ dinv_ed[i].first ] = dinv_ed[i].second;

               // Apply graph edits.

               if (verbose)
               {    cout << Date( ) << ": sort edits, peak mem = " 
                         << PeakMemUsageGBString( ) << endl;    }
               ParallelSort(edits);
               if (verbose)
               {    cout << Date( ) << ": insert edits, peak mem = " 
                         << PeakMemUsageGBString( ) << endl;    }
               for ( int64_t i = 0; i < edits.isize( ); i++ )
               {    int pv1 = edits[i].first, pv2 = edits[i].second; 
                    int pe = edits[i].third;
                    // CAN THIS HAPPEN AND WHAT DOES IT MEAN?
                    if ( Member( D.FromEdgeObj(pv1), pe ) ) continue;
                    D.FromMutable(pv1).push_back(pv2); 
                    D.ToMutable(pv2).push_back(pv1);
                    D.FromEdgeObjMutable(pv1).push_back(pe); 
                    D.ToEdgeObjMutable(pv2).push_back(pe);    }    }    }
     if (verbose) cout << Date( ) << ": sortsyncing" << endl;
     #pragma omp parallel for
     for ( int i = 0; i < nvreps; i++ )
     {    SortSync( D.FromMutable(i), D.FromEdgeObjMutable(i) );
          SortSync( D.ToMutable(i), D.ToEdgeObjMutable(i) );    }
     if (verbose)
     {    cout << Date( ) << ": initialization complete, peak mem = " 
               << PeakMemUsageGBString( ) << endl;    }

     // Check inversion.

     int64_t unset = 0, uninv = 0;
     for ( int64_t i = 0; i < dinv.jsize( ); i++ )
     {    if ( dinv[i] < 0 ) unset++;
          else if ( dinv[dinv[i]] != i ) uninv++;    }
     ForceAssertEq( unset, 0 );
     ForceAssertEq( uninv, 0 );    }

template vec<int>& digraphE<vec<int>>::EdgeObjectMutable(int);
