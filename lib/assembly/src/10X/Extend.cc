// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS


#include "CoreTools.h"
#include "10X/Extend.h"
#include "feudal/PQVec.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/GapToyTools.h"

template <class VQ>
void ExtendPathsNew( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& bases, const VQ& qualsx, ReadPathVec& paths,
     const Bool BACK_EXTEND )
{
     // Ancillary stuff.

     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     int K = hb.K( );

     // Extend path backwards if possible.

     if (BACK_EXTEND)
     {    cout << Date( ) << ": extending paths backward" << endl;
          #pragma omp parallel for schedule(dynamic, 10000)
          for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          {    ReadPath& p = paths[id];
               int offset = p.getOffset( );
               if ( p.size( ) == 0 ) continue;
               while( offset < 0 && hb.To( to_left[ p[0] ] ).nonempty( ) )
               {    int v = to_left[ p[0] ];
                    Bool extended = False;
                    for ( int j = 0; j < hb.To(v).isize( ); j++ )
                    {    int e = hb.ITo( v, j );
                         Bool mismatch = False;
                         int n = hb.Kmers(e);
                         for ( int l = 0; l < hb.Bases(e); l++ )
                         {    if ( l - offset - n < 0 ) continue;
                              if ( hb.O(e)[l] != bases[id][ l - offset - n ] )
                              {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;
                         extended = True;
                         p.push_front(e);
                         p.setOffset( p.getOffset( ) + n );
                         offset += n;
                         break;    }
                    if ( !extended ) break;    }    }    }

     // Extend paths 1 (to right).

     const int MIN_GAIN = 5;
     const int EXT_MODE = 1;
     cout << Date( ) << ": extending paths" << endl;
     auto qualsx_clone = qualsx;
     #pragma omp parallel for firstprivate(qualsx_clone)
     for ( int64_t i = 0; i < (int64_t) paths.size( ); i++ )
     {    if ( paths[i].size( ) > 0 ) paths[i].resize(1);
          ExtendPath( paths[i], i, hb, to_right, bases[i], qualsx_clone[i],
                  MIN_GAIN, False, EXT_MODE );    }
     Validate( hb, inv, paths );

     // Extend paths 2.

     cout << Date( ) << ": starting secondary path extension" << endl;
     #pragma omp parallel for
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
     {    ReadPath& p = paths[id];
          if ( p.size( ) == 0 ) continue;

          // Do last K bases of path match perfectly?

          int pp = 0;
          int e = p[pp];
          int rpos = 0, epos = p.getOffset( );
          if ( epos < 0 )
          {    rpos = -epos;
               epos = 0;    }
          int mismatch = -1;
          const basevector& b = bases[id];
          while(1)
          {    if ( b[rpos] != hb.O(e)[epos] )
               {    mismatch = rpos;
                    break;   }
               rpos++;
               if ( rpos == b.isize( ) ) break;
               epos++;
               if ( epos == hb.O(e).isize( ) )
               {    pp++;
                    if ( pp == (int) p.size( ) ) break;
                    e = p[pp];
                    epos = K - 1;    }    }
          if ( mismatch >= 0 && mismatch > rpos - K ) continue;

          // Extend.

          if ( rpos == b.isize( ) ) continue;
          while(1)
          {    if ( epos < hb.O(e).isize( ) )
               {    if ( b[rpos] != hb.O(e)[epos] ) break;    }
               else
               {    int v = to_right[e];
                    Bool ext = False;
                    for ( int j = 0; j < hb.From(v).isize( ); j++ )
                    {    int f = hb.IFrom( v, j );
                         if ( b[rpos] == hb.O(f)[K-1] )
                         {    e = f;
                              epos = K-1;
                              p.push_back(e);
                              ext = True;    }    }
                    if ( !ext ) break;    }
               rpos++;
               epos++;
               if ( rpos == b.isize( ) ) break;    }

          /*
          if ( (int) p.size( ) > np )
          {
               #pragma omp critical
               {    cout << "\nid = " << id << ", p = " << printSeq(p) << endl;
                    cout << "inv(p) = ";
                    for ( int j = (int) p.size( ) - 1; j >= 0; j-- )
                    {    if ( j < (int) p.size( ) - 1 ) cout << ",";
                         cout << inv[ p[j] ];    }
                    cout << endl;    }    }
          */
          }

     Validate( hb, inv, paths );    }

template void ExtendPathsNew( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& bases, const VecPQVec& qualsx, ReadPathVec& paths,
     const Bool BACK_EXTEND );

template void ExtendPathsNew( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& bases, const VirtualMasterVec<PQVec>& qualsx, ReadPathVec& paths,
     const Bool BACK_EXTEND );
