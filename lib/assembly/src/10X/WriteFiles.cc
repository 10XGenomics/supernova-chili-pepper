// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Intvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "10X/DfTools.h"
#include "10X/WriteFiles.h"
#include "10X/astats/GenomeAlign.h"

void WriteAssemblyFiles( const HyperBasevector& hb, const vec<int>& inv,
     ReadPathVec& paths, VecULongVec& paths_index, const vec<int64_t>& bci, 
     const Bool ALIGN, const vecbasevector& genome, const String& dir,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb )
{
     cout << Date( ) << ": writing files" << endl;
     cout << Date( ) << ": hb has checksum " << hb.CheckSum( ) << endl;
     Mkdir777(dir);
     if ( paths.size() > 0 ) {
          cout << Date() << ": writing paths" << endl;
          paths.WriteAll( dir + "/a.paths" );
          paths.resize(0);
     }
     Echo( ToString( hb.K( ) ), dir + "/a.k" );
     BinaryWriter::writeFile( dir + "/a.hbv", hb );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     BinaryWriter::writeFile( dir + "/a.to_left", to_left );
     BinaryWriter::writeFile( dir + "/a.to_right", to_right );
     BinaryWriter::writeFile( dir + "/a.inv", inv );
     HyperBasevectorX hbx(hb);
     BinaryWriter::writeFile( dir + "/a.hbx", hbx );
     cout << Date() << ": Before removing paths and paths_index, mem = " 
          << MemUsageGBString() << endl;
     {
          // Write paths index.

          if ( paths_index.size() > 0 ) {
               cout << Date( ) << ": writing paths index" << endl;
               // the way paths_index is computed by invert(), if the highest value
               // edge has no read support, then it will be excluded.  This is bad.
               ForceAssertEq( paths_index.size(), hb.E() );
               paths_index.WriteAll( dir + "/a.paths.inv" );
               paths_index.resize(0);
          }

          cout << Date() << ": Destroyed paths and paths_index, mem = " 
               << MemUsageGBString() << endl;

          // Write edges and kmers.

          cout << Date( ) << ": making and writing edges" << endl;
          {    vecbvec edges( hb.Edges( ).begin( ), hb.Edges( ).end( ) );
               edges.WriteAll( dir + "/a.fastb" );    }
          {    vec<int> kmers( hb.E( ) );
               for ( int e = 0; e < hb.E( ); e++ )
                    kmers[e] = hb.Kmers(e);
               BinaryWriter::writeFile( dir + "/a.kmers", kmers );    }

          // Align to genome.

          if (ALIGN)
          {    const int align_genome_K = 80;
               GenomeAlign<align_genome_K>( hbx, inv, genome, alignsb );
               alignsb.WriteAll( dir + "/a.alignsb" );    }
          else Remove( dir + "/a.alignsb" );    }

     cout << Date() << ": reading paths and paths_index" << endl;
     paths_index.ReadAll( dir + "/a.paths.inv");
     {
          const int ns = 1;
          vec<vec<int>> count( ns, vec<int>( hb.E( ), 0 ) );
#pragma omp parallel for
          for ( int i = 0; i < hb.E(); ++i ) {
               count[0][i] = paths_index[i].size();
               if ( inv[i] != i ) count[0][i] += paths_index[inv[i]].size();
          }
          BinaryWriter::writeFile( dir + "/a.countsb", count );
     }

     paths.ReadAll( dir + "/a.paths");
}
