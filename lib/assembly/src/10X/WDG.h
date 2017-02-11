// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_WDG_H
#define TENX_WDG_H

#include "CoreTools.h"
#include "graph/Digraph.h"

class om;

// WDG.  Given a unique-sorted collection of words w, imagine that the words define
// a digraph consisting of lines.  For example two words make two lines, like this:
//
//        --A--> --X--> --B-->
//        --A--> --Y--> --B-->
//
// Now suppose also given a set of overlaps between those words.  Then create the
// identification (quotient) graph obtained by squeezing those words together.
// For example if we overlap both copies of A, and do like wise for B, then we get
//
// four edges (poorly rendered here):
//
//        --A--> --X--> --B-->
//               --Y-->
//
// To save space, the algorithm deletes its inputs when done with them.

void WDG( vec< vec<int> >& w, vec<int64_t>& cinv, vec<vec<om>>& omatch, 
     digraphE<int>& D, vec<int>& dinv, const Bool verbose );

class om {

     public:

     om( ) { }
     om( const int c2, const int start1, const int start2, const int len ) :
          c2(c2), start1(start1), start2(start2), len(len) { }

     int c2;
     int start1, start2;
     int len;

     int Offset( ) const { return start1 - start2; }

     int Start1( ) const { return start1; }
     int Start2( ) const { return start2; }
     int Stop1( ) const { return start1 + len; }
     int Stop2( ) const { return start2 + len; }

     void Extend( const vec<int>& C1, const vec<int>& C2 )
     {    int n1 = C1.size( ), n2 = C2.size( );
          while( start1 > 0 && start2 > 0 && C1[start1-1] == C2[start2-1] )
          {    start1--;
               start2--;
               len++;    }
          while( start1 + len < n1 && start2 + len < n2
               && C1[start1+len] == C2[start2+len] )
          {    len++;    }    }

     void Validate( const vec<int>& C1, const vec<int>& C2 )
     {    for ( int i = 0; i < len; i++ )
          {    if ( C1[ Start1( ) + i ] != C2[ Start2( ) + i ] )
               {    cout << "\nInvalid om." << endl;
                    TracebackThisProcess( );    }    }    }

     friend Bool operator==( const om& o1, const om& o2 )
     {    return o1.c2 == o2.c2 && o1.start1 == o2.start1 && o1.start2 == o2.start2 
               && o1.len == o2.len;    }

};

TRIVIALLY_SERIALIZABLE(om);

#endif
