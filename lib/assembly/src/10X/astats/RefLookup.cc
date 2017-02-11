// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "FetchReads.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "random/Shuffle.h"
#include "10X/astats/RefLookup.h"

void FetchFinished( const String& SAMPLE, vecbasevector& G )
{    if ( SAMPLE == "NA12878" )
          FetchReads( G, 0, "/mnt/opt/meowmix_git/assembly/refs/fos100.fasta" );
     else if ( SAMPLE == "HGP" )
     {    String genome = "/mnt/home/jaffe/refs/mrx/genome.fastb";
          int n = MastervecFileObjectCount(genome);
          vec<int> x;
          Shuffle( n, x );
          x.resize(400);
          G.Read( genome, x );    }    }

template<int K> void RefLookup( const HyperBasevectorX& hb, 
     const vecbasevector& G, vec< vec< pair<int,int> > >& locs )
{
     // Build reference lookup table.

     cout << Date( ) << ": building lookup table for reference" << endl;
     ForceAssertEq( K, hb.K( ) );
     vec< triple<kmer<K>,int,int> > gkmers_plus;
     MakeKmerLookup0( G, gkmers_plus );

     // Find the reference kmers that are in the assembly.

     locs.clear_and_resize( G.size( ) );
     {    for ( int g = 0; g < (int) G.size( ); g++ )
               locs[g].resize( G[g].isize( ) - K + 1, make_pair( -1, -1 ) );
          const int64_t batches = 1000;
          int N = hb.E( );
          cout << Date( ) << ": finding reference kmers in the assembly" << endl;
          vec<vec<quad<int,int,int,int>>> locsb( batches);
          double fclock = WallClockTime( );
          #pragma omp parallel for schedule(dynamic, 1)     // one batch per thread
          for ( int64_t b = 0; b < batches; b++ )
          {    kmer<K> x;
               // for edges in the batch
               for ( int e = (b*N)/batches; e < ((b+1)*N)/batches; e++ )
               {    const basevector& E = hb.O(e);
                    // kmerize the edge
                    for ( int epos = 0; epos <= E.isize( ) - K; epos++ )
                    {    x.SetToSubOf( E, epos );
                         int64_t low = LowerBound1( gkmers_plus, x );
                         if ( low < gkmers_plus.jsize( ) 
                              && gkmers_plus[low].first == x )
                         {    int64_t high;
                              for ( high = low + 1; high < gkmers_plus.jsize( ); 
                                   high++ )
                              {    if ( gkmers_plus[high].first != x ) break;    }
                              // kmers [low, high) match this kmer
                              for ( int64_t j = low; j < high; j++ )
                              {    int g = gkmers_plus[j].second; 
                                   int gpos = gkmers_plus[j].third;
                                   // locs[genome contig][genome position] <-
                                   //      ( edgeno, start position on the edge )
                                   locsb[b].push_back( 
                                        make_quad( g, gpos, e, epos ) );
                                             }    }    }    }
               Sort( locsb[b] );    }
          for ( int64_t b = 0; b < batches; b++ )
          for ( int64_t j = 0; j < locsb[b].jsize( ); j++ )
          {    int g = locsb[b][j].first;
               int gpos = locsb[b][j].second;
               int e = locsb[b][j].third;
               int epos = locsb[b][j].fourth;
               locs[g][gpos] = make_pair( e, epos );    }    }    }

template void RefLookup<40>( const HyperBasevectorX& hb, 
     const vecbasevector& G, vec< vec< pair<int,int> > >& locs );

template void RefLookup<48>( const HyperBasevectorX& hb, 
     const vecbasevector& G, vec< vec< pair<int,int> > >& locs );

template void RefLookup<60>( const HyperBasevectorX& hb, 
     const vecbasevector& G, vec< vec< pair<int,int> > >& locs );
