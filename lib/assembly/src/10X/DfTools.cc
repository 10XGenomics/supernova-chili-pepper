// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Intvector.h"
#include "ParallelVecUtilities.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/ExtractReads.h"
#include "10X/DfTools.h"
#include "10X/Gap.h"

PerfStatLoggerM PerfStatLoggerM::gInst;

vec<int> GetBarcodes( const int e, const vec<int>& inv,
     const VecULongVec& paths_index, const vec<int>& bc )
{    vec<int> bs;
     for ( int x : { e, inv[e] } )
     {    for ( int j = 0; j < (int) paths_index[x].size( ); j++ )
          {    int64_t id = paths_index[x][j];
               if ( bc[id] > 0 ) bs.push_back( bc[id] );    }    }
     UniqueSort(bs);
     return bs;    }

void LoadData( const String& work_dir, const String& R, const vec<String>& lr,
     const vec<double>& LR_SELECT_FRAC, vecbasevector& bases,
     ObjectManager<VecPQVec>& quals_om, vec<int64_t>& bci,
     vec<String>& subsam_names, vec<int64_t>& subsam_starts, vec<DataSet>& datasets )
{
     String SAMPLE, species;
     String tmp_dir1 = work_dir + "/data";
     auto& quals = quals_om.object( );

		// Load ordinary read data.

		if ( R.size( ) > 0 ) {
			vec<String> regions;
			String SELECT_FRAC = "1";
			ExtractReads( SAMPLE, species, R, SELECT_FRAC, -1,
				regions, tmp_dir1, work_dir, False, False, False,
				subsam_names, subsam_starts, &bases, quals_om );
               // for now, we're assuming R is a single, PCR-free dataset
               datasets.push_back( { ReadDataType::PCR_FREE, 0 } );
		}

		bci.push_back(0);	// null barcode always at the front


		// load linked-read (LR) data in two passes to avoid duplicating
		// reads in memory -- first unbarcoded, then barcoded
          //
          if ( lr.size() ) ForceAssertEq(lr.size(), LR_SELECT_FRAC.size() );
		cout << Date() << ": reading in linked read data" << endl;
		enum passes:int { PASS_UNBARCODED, PASS_BARCODED, PASS_LAST };
		for ( int pass = PASS_UNBARCODED; pass < PASS_LAST; ++pass ) {

			for ( size_t i = 0; i < lr.size(); ++i ) {

				String head = lr[i].Before(".fastb");
				VirtualMasterVec<basevector> basesi( head + ".fastb" );
				VirtualMasterVec<PQVec> qualsi( head + ".qualp" );

				vec<int64_t> bcii;
				BinaryReader::readFile( head + ".bci" , &bcii );

				// sanity check .bci to avoid old code
				if ( bcii[0] != 0 )
					FatalErr("barcode 0 is unbarcoded data and must start at 0");

                    auto frac = LR_SELECT_FRAC[i];
				auto decider = [frac]() {
					return (1. * randomx() / RNGen::RNGEN_RAND_MAX) <= frac;
				};

				ForceAssertEq( bcii[1] % 2, 0 );
				ForceAssertLe( bcii[1], basesi.size() );
				ForceAssertLe( bcii[1], qualsi.size() );

				// TODO: these two cases below can now be merged easily
                    int64_t nbases = bases.size();
				if ( pass == PASS_UNBARCODED ) {
                         datasets.push_back( { ReadDataType::UNBAR_10X, nbases } );

					// append bases, quals from [ 0, bcii[1] )
					// assume now that we only have pairs
					for ( size_t i = 0; i < (size_t) bcii[1]; i+=2 ) {
						if ( decider() ) {
							bases.push_back(basesi[i]);
							bases.push_back(basesi[i+1]);
							quals.push_back(qualsi[i]);
							quals.push_back(qualsi[i+1]);
						}
					}

				} else if ( pass == PASS_BARCODED ) {
                         datasets.push_back( { ReadDataType::BAR_10X, nbases } );

					for ( size_t bc = 1; bc < bcii.size()-1; ++bc ) {
						// for each barcode
						bci.push_back( bases.size() );
						for ( size_t i = (size_t) bcii[bc]; i < (size_t) bcii[bc+1]; i+=2 ) {
							// for each pair of this barcode
							if ( decider() ) {
								bases.push_back(basesi[i]);
								bases.push_back(basesi[i+1]);
								quals.push_back(qualsi[i]);
								quals.push_back(qualsi[i+1]);
							}
						}
					}

				} else FatalErr("Bug - barcodes pass");

			}	// for lr...

		} // for pass...

		bci.push_back( bases.size() );

		bases.WriteAll( work_dir + "/data/frag_reads_orig.fastb" );
		quals_om.store();
		ForceAssertEq( bases.size( ), quals.size( ) );
		BinaryWriter::writeFile( work_dir + "/data/frag_reads_orig.bci", bci );

	}

void GetQualStats( const VecPQVec& quals, vec<vec<vec<int64_t>>>& hist, 
                   int & max_read_length )
{    const int64_t batch = 100000;
     int64_t N = quals.size( );
     // compute the histogram
     // hist is a 2 x max_read_length x 256 histogram
     cout << Date( ) << ": computing quality histogram" << endl; 
     // first compute the max quality score.
     /*unsigned char max_qual=0;
     #pragma omp parallel for schedule (dynamic, batch)
     for ( int64_t id = 0; id < (int64_t) quals.size( ); id ++ )
     {    qualvector q;
          quals[id].unpack(&q);
          #pragma omp critical
          {    
               for (auto x: q) {
                    if ( max_qual < x )
                         max_qual = x;
               }
          }
     }*/
     // create the data structure
     hist.resize( 2, vec<vec<int64_t>>( max_read_length, vec<int64_t>(256, 0) ) );
     #pragma omp parallel for schedule (dynamic, batch)
     for ( int64_t id = 0; id < (int64_t) quals.size( ); id ++ )
     {    qualvector q;
          quals[id].unpack(&q);
          #pragma omp critical
          {    
               int pos = 0;
               for ( auto x : q ) {
                    hist[id%2][pos][x]++;
                    pos++;
               }
          }    
     }
     // find max quality in data structure.
     int max_qual;
     bool found_max=false;
     for (max_qual = 255; max_qual >= 0; max_qual--) {
          for (int pos = 0; pos != max_read_length; pos++) {
               if ( hist[0][pos][max_qual] > 0 || hist[1][pos][max_qual] > 0 ) {
                    found_max = true;
                    break;
               }
          }
          if ( found_max )
               break;
     }
     // resize data structure
     for ( int pos = 0; pos != max_read_length; pos++ ) {
          hist[0][pos].resize(max_qual+1);
          hist[1][pos].resize(max_qual+1);
     }
     
}

void FragDist( const HyperBasevectorX& hb, const vec<int>& inv,
     const ReadPathVec& paths, vec<int64_t>& count )
{    const int max_sep = 1000;
     count.resize( max_sep + 1, 0 );
     for ( int64_t id1 = 0; id1 < (int64_t) paths.size( ); id1 += 2 )
     {    int64_t id2 = id1 + 1;
          if ( paths[id1].size( ) == 0 || paths[id2].size( ) == 0 ) continue;
          int e1 = paths[id1][0], e2 = inv[ paths[id2][0] ];
          int epos1 = paths[id1].getOffset( );
          if ( e1 != e2 ) continue;
          int n = hb.Bases(e1);
          if ( epos1 + max_sep > n ) continue;
          int epos2 = hb.Bases(e2) - paths[id2].getOffset( );
          int len = epos2 - epos1;
          if ( len < 0 || len > max_sep ) continue;
          count[len]++;    }    }

void ReadTwoPctProper( const HyperBasevectorX& hb, const vec<int>& inv,
     VirtualMasterVec<ReadPath> const& vpaths, double & r2_pct_proper )
{
     cout << Date( ) << ": computing read two pct proper";
     const int max_sep = 1000;
     double r2_map = 0.0;
     int64_t sample = 0;
     auto paths = vpaths;
     #pragma omp parallel for reduction(+:r2_map,sample) firstprivate(paths)
     for ( int64_t id1 = 0; id1 < (int64_t) paths.size( ); id1 += 2 ) {
          int64_t id2 = id1 + 1;
          // new stuff here
          // R1 must be mapped to edge e1
          if ( paths[id1].size( ) == 0 ) continue;
          int e1 = paths[id1][0];
          int le1 = hb.Bases(e1);
          int epos1 = paths[id1].getOffset( );
          
          // is the edge long enough for R2 to map here
          if ( le1 - epos1 >= max_sep ) {
               // did R2 map?
               bool mapped = false;
               if ( paths[id2].size( ) > 0 ) 
               {
                    int e2 = inv[ paths[id2][0] ];
                    if ( e1 == e2 ) {
                         //R1,R2 mapped to the same edge
                         int epos2 = hb.Bases(e2) - paths[id2].getOffset( );
                         int len = epos2 - epos1;
                         // make sure the insert isn't too long
                         if ( len >= 0 && len <= max_sep )
                              mapped=true;
                    }
               }
               if (mapped)
                    r2_map++;
               sample++;
          }
     }
     cout << " [sample size " << sample << "]" << endl;
     r2_pct_proper = (r2_map/sample)*100; // This number is a percentage!

}
template<int K> void MapClosures(
     const HyperBasevectorX& hb, const vec<basevector>& closures, ReadPathVec& clop )
{    cout << Date( ) << ": building closure lookup" << endl;

     // Remove short closures.  This should have been done elsewhere!

     vec<Bool> to_delete( closures.size( ), False );
     vec<basevector> closures2(closures);
     for ( int i = 0; i < closures.isize( ); i++ )
          if ( closures[i].isize( ) < K ) to_delete[i] = True;
     EraseIf( closures2, to_delete );

     clop.resize( closures2.size( ) );
     vec< pair<kmer<K>,int> > X( closures2.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < closures2.isize( ); i++ )
     {    const basevector& c = closures2[i];
          kmer<K> x;
          x.SetToSubOf( c, 0 );
          X[i] = make_pair( x, i );    }
     ParallelSort(X);
     cout << Date( ) << ": mapping closures back to assembly" << endl;
     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int g = 0; g < hb.E( ); g++ )
     {    const basevector& E = hb.O(g);
          kmer<K> x;
          for ( int r = 0; r < hb.Kmers(g); r++ )
          {    x.SetToSubOf( E, r );
               int64_t low = LowerBound1( X, x ), high = UpperBound1( X, x );
               for ( int64_t m = low; m < high; m++ )
               {    int c = X[m].second;
                    const basevector& C = closures2[c];
                    int e = g, epos = r;

                    // closure c starts at position epos on on edge e
                    // trace it through the assembly

                    int offset = epos;
                    vec<int> p = {e};
                    int cpos = 0;
                    while(1)
                    {    cpos += hb.Bases(e) - epos;
                         if ( cpos >= C.isize( ) ) break;
                         int v = hb.ToRight(e);
                         Bool extended = False;
                         for ( int j = 0; j < (int) hb.From(v).size( ); j++ )
                         {    int f = hb.IFrom( v, j );
                              const basevector& F = hb.O(f);
                              if ( F[K-1] == C[cpos] )
                              {    extended = True;
                                   p.push_back(f);
                                   e = f;
                                   epos = K-1;
                                   break;    }    }
                         if ( !extended ) break;    }
                    ForceAssertGe( cpos, C.isize( ) );
                    clop[c].setOffset(offset);
                    for ( int j = 0; j < p.isize( ); j++ )
                         clop[c].push_back( p[j] );    }    }    }    }

template void MapClosures<40>( const HyperBasevectorX&, const vec<basevector>&,
     ReadPathVec& );
template void MapClosures<48>( const HyperBasevectorX&, const vec<basevector>&,
     ReadPathVec& );
template void MapClosures<60>( const HyperBasevectorX&, const vec<basevector>&,
     ReadPathVec& );

void SanityCheckBarcodeCounts( const vec<int64_t>& bci )
{    if ( bci.size( ) < 2 ) return;
     const int big = 50000;
     const double max_big_frac = 0.1;
     int64_t big_total = 0;
     for ( int j = 2; j < bci.isize( ); j++ )
     {    int64_t n = bci[j] - bci[j-1];
          if ( n >= big ) big_total += n;    }
     int64_t total = bci.back( ) - bci[1];
     if ( total > 0 && double(big_total) / double(total) > max_big_frac )
     {    cout << "\nOf your barcoded reads, "
               << PERCENT_RATIO( 3, big_total, total ) << " lie in "
               << "barcodes having at least " << big << " reads." << endl;
          cout << "That would suggest that something is very wrong." << endl;
          cout << "Giving up.\n" << endl;
          Scram(1);    }    }

void MakeDots( int& done, int& ndots, const int total )
{    if ( done % ( (total+99) / 100 ) == 0 )
     {    if ( ndots < 100 )
          {    cout << ".";
               ndots++;
               if ( ndots > 0 && ndots % 50 == 0 ) cout << "\n";
               else if ( ndots > 0 && ndots % 10 == 0 ) cout << " ";
               flush(cout);    }    }
     if ( done == total - 1 )
     {    while ( ndots < 100 )
          {    cout << ".";
               ndots++;
               if ( ndots > 0 && ndots % 50 == 0 ) cout << "\n";
               else if ( ndots > 0 && ndots % 10 == 0 ) cout << " ";
               flush(cout);    }    }
     done++;    }

// Note grossly inefficient conversion below.
// Note that we might want to factor this through Munch.

void SuperToSeqGraph( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     HyperBasevectorX& hbd )
{    vec<basevector> edges( D.E( ) );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<int> seqverts;

     // Convert most edges.

     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] >= 0 ) edges[d] = hb.Cat( D.O(d) );
          else if ( IsSequence( D.O(d) ) ) 
          {
               #pragma omp critical
               {    seqverts.push_back( to_left[d] );    }    }    }
     UniqueSort(seqverts);

     // Go through the vertices to left of sequence gaps.
     // Not necessarily symmetric.

     int nprobs = 0;
     for ( int i = 0; i < seqverts.isize( ); i++ )
     {    int v = seqverts[i];
          int ltrim, rtrim;

          // Test for problem.

          int d = D.IFrom( v, 0 );
          int w = to_right[d];
          GapToSeq( D.O(d), ltrim, rtrim, edges[d] );    
          if ( D.To(v).empty( ) || D.From(w).empty( ) )
          {    cout << "\nPROBLEM:\n";
               PRINT3( d, D.To(v).size( ), D.From(w).size( ) );
               Scram(0);    }
          int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
          Bool problem = False;
          if ( edges[d1].isize( ) - ltrim < hb.K( ) )
          {    // cout << "\nPROBLEM:\n";
               // PRINT4( d, d1, ltrim, edges[d1].size( ) );
               problem = True;    }
          if ( edges[d2].isize( ) - rtrim < hb.K( ) )
          {    // cout << "\nPROBLEM:\n";
               // PRINT4( d, d2, rtrim, edges[d2].size( ) );
               problem = True;    }
          if ( d1 == d2 && edges[d1].isize( ) - ltrim - rtrim < hb.K( ) )
               problem = True;
          if (problem)
          {    for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    int d = D.IFrom( v, j );
                    edges[d].clear( );    }
               nprobs++;
               continue;    }

          // Proceed.

          for ( int j = 0; j < D.From(v).isize( ); j++ )
          {    int d = D.IFrom( v, j );
               int w = to_right[d];
               if ( !D.To(v).solo( ) ) cout << "Problem 1 at edge " << d << endl;
               ForceAssert( D.To(v).solo( ) );
               if ( !D.From(w).solo( ) ) cout << "Problem 2 at edge " << d << endl;
               ForceAssert( D.From(w).solo( ) );
               GapToSeq( D.O(d), ltrim, rtrim, edges[d] );    
               if ( j == 0 )
               {    int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
                    if ( ltrim > 0 ) edges[d1].resize( edges[d1].isize( ) - ltrim );
                    if ( rtrim > 0 ) 
                    {    if ( edges[d2].isize( ) - rtrim < hb.K( ) )
                         {    cout << "\nPROBLEM:\n";
                              PRINT4( d, d2, rtrim, edges[d2].size( ) );    }
                         edges[d2].SetToSubOf( edges[d2], rtrim, 
                              edges[d2].isize( ) - rtrim );    }    }    }    }
#ifndef CS
     cout << Date( ) << ": problems converting " << nprobs << " of "
          << seqverts.size( ) << " gaps" << endl;
#endif
     HyperBasevector H( hb.K( ), D, edges );
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     hbd = HyperBasevectorX(H);    }
