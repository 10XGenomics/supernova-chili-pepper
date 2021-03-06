// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Note that much of the code here that assays assemblies using a reference
// sequence is designed for HUMAN samples.  All of it would have to be carefully
// reviewed and modified to allow for nonhuman samples.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "10X/DfTools.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/Super.h"
#include "10X/astats/LineLine.h"
#include "10X/astats/MeasureGaps.h"
#include "10X/astats/Misassembly.h"
#include "10X/astats/N50P.h"
#include "10X/astats/View.h"

template<class T>
void EraseIf( SerfVec<T>& v, const vec<Bool>& to_delete )
{    SerfVec<T> v2;
     for ( int64_t i = 0; i < (int64_t) v.size( ); i++ )
          if ( !to_delete[i] ) v2.push_back( v[i] );
     v = v2;    }

void ReportAssemblyStats( const vecbasevector& genome,
     const vec< pair<int,ho_interval> >& ambint, const HyperBasevectorX& hb,
     const vec<int>& inv, digraphE<vec<int>> D, vec<int> dinv,
     const vec<vec<vec<vec<int>>>>& dlines, 
     const ReadPathVec& dpaths, const vec<double>& COV, 
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb, const vecbasevector& G, 
     const vec<vec<pair<int,int>>>& locs, double & r2_pct_proper, ostream& out, 
     const String& DIR, const String& WRITE_SUB,
     const String& CS_SAMPLE_ID, const String& CS_SAMPLE_DESC )
{
     // Get checksum.

     auto checksum = CheckSum( hb, inv, D, dinv );

     // Declare the minimum line that we'll use.

     const int MIN_LINE = 10000;

     // Compute some auxiliary stuff.

     vec<int> kmers( hb.E( ) ), llens;
     // We have two versions of GetLineLengths, one with order D, hb
     // and one with order hb, D.  I hope these are identical.  The order
     // used below also appears in one other place in this file.
     GetLineLengths( D, hb, dlines, llens );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          kmers[e] = hb.Kmers(e);
     int K = hb.K( );

     // Get linelocs.

     vec<vec< pair<int,int> >> linelocs( kmers.size( ) );
     for ( int i = 0; i < dlines.isize( ); i++ )
     for ( int j = 0; j < dlines[i].isize( ); j++ )
     for ( int k = 0; k < dlines[i][j].isize( ); k++ )
     for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
     {    int d = dlines[i][j][k][l];
          if ( D.O(d)[0] < 0 ) continue;
          for ( int m = 0; m < D.O(d).isize( ); m++ )
          {    int e = D.O(d)[m];
               linelocs[e].push( i, j );    }    }

     // Track misassembly stats.

     int64_t total_err_dis_num = 0, total_err_dis_den = 0;
     int64_t total_err_ori_num = 0, total_err_ori_den = 0;
     int64_t total_err_ord_num = 0, total_err_ord_den = 0;

     // Compute coverage.  Horrible, done by deleting "one chromosome".

     auto per = [&]( double n, double d )
     {    ostringstream out;
          out << fixed << setprecision(2) << setw(7) << right << 100*n/d;
          return out.str( );    };
     int64_t genome_size = 0, cap_gap_total = 0, cov_total = 0;
     int64_t contig_line_N50 = 0, pairtig_N50 = 0;
     String OUTDIR = DIR + "/" + "final" + WRITE_SUB;

     // Find lines of lines.

     cout << Date( ) << ": finding lines of lines" << endl;
     vec<vec<vec<vec<int>>>> dlines2;
     FindLineLines( D, dinv, dlines, dlines2 );
     BinaryWriter::writeFile( OUTDIR + "/a.sup.linelines", dlines2 );
     vec<int> linv;
     LineInv( dlines, dinv, linv );

     // For each line line bubble, delete one branch.

     vec<int> dels;
     for ( int i = 0; i < dlines2.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines2[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    if ( L[j].size( ) == 2 && L[j][0].solo( ) && L[j][1].solo( ) )
               {    vec<int> x = Contents( L[j] );
                    if ( x.solo( ) ) continue; // should never happen
                    vec<int> y = { linv[ L[j][0][0] ], linv[ L[j][1][0] ] };
                    Sort(y);
                    if ( y < x ) continue;
                    if ( llens[ x[0] ] >= llens[ x[1] ] )
                         dels.append( Contents( dlines[ x[1] ] ) );
                    else dels.append( Contents( dlines[ x[0] ] ) );    }    }    }
     int nd = dels.size( );
     for ( int i = 0; i < nd; i++ ) dels.push_back( dinv[ dels[i] ] );
     digraphE<vec<int>> F(D);
     vec<int> finv(dinv);
     cout << Date( ) << ": deleting edges from F" << endl;
     F.DeleteEdgesParallel(dels);
     RemoveUnneededVertices( F, finv );
     CleanupCore( F, finv );

     // Find lines, and locations of base edges on them.

     vec<vec<vec<vec<int>>>> flines;
     cout << Date( ) << ": finding lines again" << endl;
     FindLines( F, finv, flines, MAX_CELL_PATHS, MAX_CELL_DEPTH, False );
     cout << Date( ) << ": computing linelocs" << endl;
     vec<vec< pair<int,int> >> linelocs2( kmers.size( ) );
     for ( int i = 0; i < flines.isize( ); i++ )
     for ( int j = 0; j < flines[i].isize( ); j++ )
     for ( int k = 0; k < flines[i][j].isize( ); k++ )
     for ( int l = 0; l < flines[i][j][k].isize( ); l++ )
     {    int d = flines[i][j][k][l];
          if ( F.O(d)[0] < 0 ) continue;
          for ( int m = 0; m < F.O(d).isize( ); m++ )
          {    int e = F.O(d)[m];
               linelocs[e].push( i, j );    }    }
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          UniqueSort( linelocs[e] );
     vec<int> llens2, to_left2, to_right2;
     GetLineLengths( F, hb, flines, llens2 );
     F.ToLeft(to_left2), F.ToRight(to_right2);

     // Compute contig N50.

     cout << Date( ) << ": computing contig N50" << endl;
     for ( int pass = 1; pass <= 2; pass++ )
     {    vec<int> lens;
          for ( int i = 0; i < flines.isize( ); i++ )
          {    if ( llens2[i] < MIN_LINE ) continue;
               int pos = 0;
               const vec<vec<vec<int>>>& L = flines[i];
               for ( int j = 0; j < L.isize( ); j++ )
               {    const vec<vec<int>>& M = L[j];

                    // Break contig at cell if every path in the cell
                    // contains either a gap (pair or barcode), or a cycle.

                    Bool gap = True;
                    for ( int k = 0; k < M.isize( ); k++ )
                    {    if ( M[k].empty( ) )
                         {    int g = -1;
                              if ( j > 0 && L[j-1].nonempty( )
                                   && L[j-1][0].nonempty( ) )
                              {    int d = L[j-1][0][0];
                                   int v = to_right2[d];
                                   if ( v >= 0 ) g = F.IFrom(v,0);    }
                              if ( IsSequence( F.O(g) ) ) gap = False;
                              if ( pass == 2 )
                              {    if ( g >= 0 && IsPairGap( F.O(g) ) )
                                        gap = False;    }    }
                         else
                         {    Bool gappy = False;
                              const vec<int>& x = M[k];
                              for ( int l = 0; l < x.isize( ); l++ )
                              {    if ( pass == 2 && IsPairGap( F.O(x[l]) ) )
                                        continue;
                                   if ( IsSequence( F.O(x[l]) ) ) continue;
                                   if ( F.O(x[l])[0] < 0 ) gappy = True;    }
                              if ( !gappy ) gap = False;    }    }
                    if (gap)
                    {    if ( pos >= 1 ) lens.push_back(pos);
                         pos = 0;    }
                    vec<int> lensj;
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    {    int len = 0;
                         for ( int l = 0; l < L[j][k].isize( ); l++ )
                         {    int d = L[j][k][l];
                              if ( IsSequence( F.O(d) ) )
                              {    int ltrim, rtrim;
                                   basevector x;
                                   GapToSeq( F.O(d), ltrim, rtrim, x );
                                   len += x.isize( ) - K + 1;
                                   int v = to_left2[d];
                                   if ( F.IFrom(v,0) == d )
                                        len = len - ltrim - rtrim;
                                   continue;    }
                              else if ( F.O(d)[0] < 0 ) continue;
                              for ( int m = 0; m < F.O(d).isize( ); m++ )
                                   len += hb.Kmers( F.O(d)[m] );    }
                         lensj.push_back(len);    }
                    Sort(lensj);
                    if ( lensj.nonempty( ) ) pos += Median(lensj);    }
               if ( pos >= 1 ) lens.push_back(pos);    }
          if ( lens.nonempty( ) )
          {    if ( pass == 1 ) contig_line_N50 = N50(lens);
               else pairtig_N50 = N50(lens);    }    }

     // Proceed with flines analysis.

     if ( genome.size( ) > 0 )
     {    
          // Quantify misassemblies.

          // #pragma omp parallel for
          for ( int L = 0; L < flines.isize( ); L++ )
          {    if ( llens2[L] < MIN_LINE ) continue;
               SerfVec< quad<int,Bool,ho_interval,ho_interval> > view;
               View( L, K, kmers, inv, F, flines, linelocs2, alignsb, view );
               Misassembly( view, total_err_dis_num, total_err_dis_den,
                    total_err_ori_num, total_err_ori_den, total_err_ord_num,
                    total_err_ord_den );    }

          // Compute coverage.

          vec< quad< int, ho_interval, int, vec<String> > > all_notes;
          cout << Date( ) << ": start main loop" << endl;
          #pragma omp parallel for schedule( dynamic, 100 )
          for ( int L = 0; L < flines.isize( ); L++ )
          {    if ( llens2[L] < MIN_LINE ) continue;
               SerfVec< quad<int,Bool,ho_interval,ho_interval> > view;
               View( L, K, kmers, inv, F, flines, linelocs2, alignsb, view );

               // Delete short and rc segments.

               vec<Bool> to_delete( view.size( ), False );
               for ( int i = 0; i < (int) view.size( ); i++ )
               {    if ( !view[i].second ) to_delete[i] = True;
                    else
                    {    int gstart = view[i].third.Start( );
                         int gstop = view[i].third.Stop( );
                         if ( gstop - gstart < MIN_LINE )
                         {    to_delete[i] = True;    }    }    }
               EraseIf( view, to_delete );

               // Merge proximate segments.

               Sort(view);
               to_delete.resize_and_set( view.size( ), False );
               for ( int i = 0; i < (int) view.size( ); i++ )
               {    int j;
                    for ( j = i + 1; j < (int) view.size( ); j++ )
                    {    if ( view[j].first != view[i].first ) break;
                         if ( view[j].third.Start( ) - view[j-1].third.Stop( )
                              > 50000 )
                         {    break;    }    }
                    view[i].third.SetStop( view[j-1].third.Stop( ) );
                    view[i].fourth.SetStop( view[j-1].fourth.Stop( ) );
                    for ( int k = i + 1; k < j; k++ ) to_delete[k] = True;
                    i = j - 1;    }
               EraseIf( view, to_delete );

               // Convert to notes.

               vec< quad< int, ho_interval, int, vec<String> > > notes;
               for ( int i = 0; i < (int) view.size( ); i++ )
               {    int gstart = view[i].third.Start( );
                    int gstop = view[i].third.Stop( );
                    int g = view[i].first;
                    int astart = view[i].fourth.Start( );
                    int astop = view[i].fourth.Stop( );
                    String chr = ToString(g+1);
                    if ( g == 22 ) chr = "X";
                    if ( g == 23 ) chr = "Y";
                    ostringstream out0, out0b, out1, out2, out3, out4, out5;
                    out0 << "L" << L;
                    // out0b << ToString( COV[L], 2 );
                    out1 << chr << ":" << ToString(gstart/1000000.0,3) + "-"
                         + ToString(gstop/1000000.0,3);
                    out2 << ToStringAddCommas(gstop-gstart);
                    out3 << ToStringAddCommas(astop-astart);

                    // Compute "gap", the estimated sum of captured gaps within
                    // the unit.

                    int64_t gap = (gstop-gstart) - (astop-astart);
                    out4 << ToStringAddCommas(gap);
                    out5 << fixed << setprecision(2) << right
                         << (100.0 * gap) / (gstop-gstart) << "%";
                    notes.push( g, view[i].third, gap, vec<String> { out0.str( ),
                         // out0b.str( ), 
                         out1.str( ), out2.str( ), out3.str( ),
                         out4.str( ), out5.str( ) } );    }
               #pragma omp critical
               {    all_notes.append(notes);    }    }
          cout << Date( ) << ": sorting notes" << endl;
          Sort(all_notes);

          // Determine if chromosome Y appears to be present.  Note different
          // treatment of ambiguous bases in the reference, because there are so
          // many of them on Y.  Data for cov_Y on which Y threshold was based:
          // NA12878 (female):  0.33%
          // NA24385 (male):   44.29%
          // HGP (male):       42.91%.
          // All these are from 1200M reads.  We set a low threshold because 
          // for lower coverage, or lower quality reads, coverage for a male sample
          // might be much lower.

          int64_t genome_size_Y_not_N = genome[23].size( );
          for ( int l = 0; l < ambint.isize( ); l++ )
          {    if ( ambint[l].first == 23 ) 
                    genome_size_Y_not_N -= ambint[l].second.Length( );    }
          int64_t cap_gap_total_Y = 0, cov_total_Y = 0;
          for ( int i = 0; i < all_notes.isize( ); i++ )
          {    int j, g = all_notes[i].first;
               if ( g != 23 ) continue;
               vec<ho_interval> cov;
               for ( j = i + 1; j < all_notes.isize( ); j++ )
                    if ( all_notes[j].first != g ) break;
               for ( int k = i; k < j; k++ )
               {    cap_gap_total_Y += Max( 0, all_notes[k].third );
                    cov.push_back( all_notes[k].second );    }
               cov_total_Y += TotalCovered(cov);
               i = j - 1;    }
          double cap_frac_Y = cap_gap_total_Y / double(genome_size_Y_not_N);
          double uncap_frac_Y 
               = (genome_size_Y_not_N-cov_total_Y) / double(genome_size_Y_not_N);
          double cov_Y = 1 - cap_frac_Y - uncap_frac_Y;
          cout << Date( ) << ": nominal coverage of chromosome Y = "
               << per( cov_Y*10000, 10000 ) << "%" << endl;
          Bool male = ( cov_Y >= 0.05 );
          cout << Date( ) << ": calling this sample "
               << ( male ? "MALE" : "FEMALE" ) << endl;
          int max_chr = ( male ? 23 : 22 );

          // Compute coverage.

          for ( int g = 0; g <= max_chr; g++ )
               genome_size += genome[g].size( );
          for ( int i = 0; i < all_notes.isize( ); i++ )
          {    int j, g = all_notes[i].first;
               if ( g > max_chr ) break;
               vec<ho_interval> cov;
               for ( j = i + 1; j < all_notes.isize( ); j++ )
                    if ( all_notes[j].first != g ) break;
               for ( int k = i; k < j; k++ )
               {    cap_gap_total += Max( 0, all_notes[k].third );
                    cov.push_back( all_notes[k].second );    }
               for ( int l = 0; l < ambint.isize( ); l++ )
                    if ( ambint[l].first == g ) cov.push_back( ambint[l].second );
               cov_total += TotalCovered(cov);
               i = j - 1;    }

          // Compute genome view.

          {    vec<vec<String>> rows;
               vec<String> row = { "xline", /* "CN", */ "genome", "glength", 
                    "alength", "gap", "gapfrac" };
               rows.push_back(row);
               rows.push_back( vec<String>( ) );
               for ( int i = 0; i < all_notes.isize( ); i++ )
               {    if ( i > 0 && all_notes[i].first != all_notes[i-1].first )
                         rows.push_back( vec<String>( ) );
                    rows.push_back( all_notes[i].fourth );    }
               Ofstream( xout, OUTDIR + "/o.genome" );
               xout << "\nWARNING: line numbers below don't make sense\n\n";
               PrintTabular( xout, rows, 3, "lrlrrrr" );    }    }

     // Measure gaps.  Partially outdated....

     double cap_frac = -1;
     if ( genome.size( ) > 0 )
     {    cout << Date( ) << ": set up to measure gaps, saving to o.gaps" << endl;
          int64_t total_lens = 0;
          int total_gaps = 0;
          for ( int i = 0; i < dlines.isize( ); i++ )
          {    total_lens += llens[i] + (K-1);
               const vec<vec<vec<int>>>& L = dlines[i];
               for ( int j = 0; j < L.isize( ); j++ )
               {    if ( !L[j].solo( ) || !L[j][0].empty( ) ) continue;
                    total_gaps++;    }    }
          vec< pair<int,int> > gaps;
          Bool verbose = False;
          MeasureGaps( K, kmers, inv, D, dlines, llens, alignsb, gaps, verbose );
          int64_t gap_sum = 0;
          for ( int i = 0; i < gaps.isize( ); i++ )
               if ( gaps[i].first >= 0 ) gap_sum += gaps[i].first;
          double total_gap_lens
               = gap_sum * double(total_gaps) / double( gaps.size( ) );
          cap_frac = total_gap_lens / total_lens;
          Ofstream( gout, OUTDIR + "/o.gaps" );
          gout << "#  gap_size  gap_edge_id";
          for ( int i = 0; i < gaps.isize( ); i++ )
               gout << gaps[i].first << "  " << gaps[i].second << "\n";    }

     // Compute total edges.

     int64_t total_edges = 0;
     for ( int e = 0; e < D.E( ); e++ )
          total_edges += D.O(e).size( );

     // Find lines of lines, and their lengths, then use this to compute 
     // scaffold N50.  Also computed estimated genome size.

     cout << Date( ) << ": finding lines of lines" << endl;
     vec<vec<vec<vec<int>>>> dlines2x;
     FindLineLines( D, dinv, dlines, dlines2x );
     vec<int> lens2;
     GetLineLineLengths( llens, dlines2x, lens2 );
     ReverseSortSync( lens2, dlines2x );
     {    Ofstream( out, OUTDIR + "/o.linelines" );
          vec<int> x;
          for ( int i = 0; i < dlines2x.isize( ); i++ )
          {    x.clear( );
               const vec<vec<vec<int>>>& M = dlines2x[i];
               for ( int j = 0; j < M.isize( ); j += 2 ) // ONLY PRINTING EVEN!
                    x.push_back( M[j][0][0] );
               out << "M" << i << "[l=" << ToStringAddCommas( lens2[i] )
                    << "]: " << printSeq(x) << "\n";    }    }
     int64_t est_genome_size = 0;
     vec<int> lllens;
     for ( auto x : lens2 ) 
     {    if ( x >= MIN_LINE ) lllens.push_back(x);
          if ( x >= MIN_LINE ) est_genome_size += x + hb.K( ) - 1;    }
     est_genome_size /= 2;
     int nlines = lllens.size( );
     int64_t scaffold_N50 = N50PL(lllens);
     int64_t scaffold_N60 = NPL(lllens, 0.6);
     cout << Date( ) << ": scaffold N50 = " << ToStringAddCommas(scaffold_N50)
          << endl;
     cout << Date( ) << ": scaffold N60 = " << ToStringAddCommas(scaffold_N60) 
          << endl;
     cout << Date( ) << ": scaffold length-weighted mean length = "
          << ToStringAddCommas( int( round( WeightedMean(lllens) ) ) ) << endl;
     cout << Date( ) << ": assembly size (sum of scaffolds >= 10 kb) = "
          << ToStringAddCommas(est_genome_size) << endl;

     // Print large scaffold sizes.  Stupid because it's showing the value for each
     // lineline and its involution, so you get pairwise duplications in the file.

     {    Ofstream( out, OUTDIR + "/stats/large_scaffold_sizes.txt" );
          for ( int i = 0; i < nlines; i++ )
               out << ToStringAddCommas( lllens[i] ) << endl;    }

     // Compute phasetig N50.

     cout << Date( ) << ": computing phasetig N50" << endl;
     vec<int> bblens; // bubble branch lengths
     vec<Bool> marked( D.E( ), False );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int i = 0; i < dlines2x.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines2x[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    if ( L[j].size( ) != 2 ) continue;
               vec<int> blens( 2, 0 );
               for ( int k = 0; k < 2; k++ )
               {    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int M = L[j][k][l];
                         blens[k] += llens[M];
                         for ( auto d : Contents(dlines[M]) )
                              marked[d] = True;    }    }
               bblens.push_back( Mean(blens) );    }    }
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    if ( L[j].size( ) < 2 ) continue;
               int d1 = L[j][0].front( ), d2 = L[j][0].back( );
               if ( marked[d1] || marked[d2] ) continue;
               int v = to_left[d1], w = to_right[d2];
               if ( D.From(v).size( ) != 2 || D.To(w).size( ) != 2 ) continue;
               int f1 = D.IFrom(v,0), f2 = D.IFrom(v,1);
               vec<int> G1, G2;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    if ( L[j][k].front( ) == f1 ) G1.push_back( L[j][k].back( ) );
                    if ( L[j][k].front( ) == f2 ) G2.push_back( L[j][k].back( ) );
                         }
               if ( !Contents(G1).solo( ) || !Contents(G2).solo( ) ) continue;
               vec<int> blens;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                         len += hb.Kmers( L[j][k][l] );
                    blens.push_back(len);    }
               Sort(blens);
               bblens.push_back( Median(blens) );    }    }
     int64_t phasetig_N50 = ( bblens.size( ) > 0 ? N50(bblens) : 0 );

     // Reinsert loops, then compute N50 perfect stretch.

     int64_t N50_perf = -1;
     double errw = -1;
     {    int nd = D.E( );
          ReinsertLoops( hb, inv, D, dinv );
          vec<Bool> keep( D.E( ), False );
          for ( int i = 0; i < dlines.isize( ); i++ )
          {    if ( llens[i] < MIN_LINE ) continue;
               vec<int> C = Contents( dlines[i] );
               for ( auto d : C ) keep[d] = True;    }
          for ( int d = nd; d < D.E( ); d++ ) keep[d] = True;
          if ( G.size( ) > 0 )
          {    String report, details;
               N50_perf = N50P( hb, D, dinv, G, locs, keep, report, details, errw );
               Ofstream( rout, OUTDIR + "/o.report" );
               rout << report;
               Ofstream( dout, OUTDIR + "/o.report.details" );
               dout << details;    }    }

     // Find N50 edge size.

     vec<int> lens;
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          int len = 0;
          for ( int i = 0; i < D.O(d).isize( ); i++ )
               len += kmers[ D.O(d)[i] ];
          lens.push_back(len);    }
     ParallelSort(lens);
     int N50_edge = N50(lens);

     // Define some printing functions.

     auto nper = []( double n, double d ) { return 100*n/d; };
     auto rat = [&]( double n, double d )
     {    ostringstream out;
          out << fixed << setprecision(2) << setw(7) << right << n/d;
          return out.str( );    };
     auto hline = [&]( )
     {    out << "----------------------------------------"
               << "----------------------------------------\n";    };
     // {    out << "========================================"
     //           << "====================================\n";    };

     // ============================================================================
     // Print final report.
     // ============================================================================

     PerfStatLoggerM::init();
     char cwd[10000];
     getcwd(cwd,10000);
     String cwds(cwd);
     cwds = cwds.RevBefore( "/" );
     cwds = cwds.RevAfter( "/" );
     out << "\n";
     hline( );
     out << "SUMMARY" << endl;
     hline( );
     out << "- " << Date( ) << endl;
     if ( CS_SAMPLE_ID == "" ) {
          PerfStatLoggerM::log( "sample_id", cwds, "pipeline sample id", true);
          out << "- [" << cwds << "]";
     } else {
          PerfStatLoggerM::log( "sample_id", CS_SAMPLE_ID, "pipeline sample id", true);
          out << "- [" << CS_SAMPLE_ID << "]";
     }
     if ( CS_SAMPLE_DESC != "" ) out << "  " << CS_SAMPLE_DESC << endl;
     if ( WRITE_SUB != "" ) out << "- SUB = " << WRITE_SUB.After( ":" ) << endl;
     // Print the code commit hash
     if ( IsRegularFile( DIR+"/../.commit_hash" ) ) {
          ifstream in_file ( DIR+"/../.commit_hash" );
          String commit_hash="";
          in_file >> commit_hash;
          PerfStatLoggerM::log( "commit_hash", commit_hash, 
               "Supernova git commit hash" );
          out << "- commit hash = " << commit_hash << endl;
     }
     PerfStatLoggerM::log( "checksum", checksum, "checksum for assembly", true);
     out << "- assembly checksum = " << ToStringAddCommas(checksum) << endl;

     // ============================================================================
     // Print input stats.
     // ============================================================================

     hline( );
     out << "INPUT" << endl;

     // Print total reads.

     int64_t nreads 
          = MastervecFileObjectCount( DIR + "/../data/frag_reads_orig.fastb" );
     out << "- " << rat(nreads,1000000) << " M   = READS"
          << "          = number of reads; ideal 800-1200 for human" <<  endl;
     PerfStatLoggerM::log( "nreads", nreads, "number of reads", true);

     // Print read length.

     {    vec<int16_t> lens;
          BinaryReader::readFile( DIR + "/../data/frag_reads_orig.lens", &lens );
          int64_t nbases = 0;
          for ( auto x : lens ) nbases += x;
          int bases_per_read = nbases / lens.size( );
          out << "- " << rat(bases_per_read,1) << " b   = MEAN READ LEN"
               << "  = mean read length after trimming; ideal 140" <<  endl;
          PerfStatLoggerM::log( 
               "bases_per_read", bases_per_read, "mean bases per read", true);    }

     // Print qual score stat.

     // test line temporary:
     if ( IsRegularFile( DIR + "/../data/frag_reads_orig.qhist" ) )
     {    vec<vec<vec< int64_t >>> qhist;
          BinaryReader::readFile( 
               DIR + "/../data/frag_reads_orig.qhist", &qhist );
          double total = 0, total30 = 0;
          int max_read_length = qhist[0].size();
          int max_qual        = qhist[0][0].size();
          for ( int pos = 0; pos < max_read_length; pos++ ) {
               for ( int qv = 0; qv != max_qual; qv++ ) {    
                    total += qhist[1][pos][qv];
                    if ( qv >= 30 ) total30 += qhist[1][pos][qv];    }
          }
          double per_Q30_R2 = 100.0 * double(total30) / double(total);
          out << "- " << rat(per_Q30_R2,1) << " %   = READ TWO Q30"
               << "   = fraction of Q30 bases in read 2; ideal 75-85" <<  endl;
          PerfStatLoggerM::log( 
               "per_q30_r2", per_Q30_R2, "read two % Q30", true );    }

     // Print insert size.

     vec<int64_t> icount;
     BinaryReader::readFile( DIR + "/a.ins_dist", &icount );
     int64_t NI = Sum(icount), isum = 0;
     int m;
     for ( m = 0; m < icount.isize( ); m++ )
     {    isum += icount[m];
          if ( isum >= NI/2 ) break;    }
     out << "- " << rat(m,1000) << " kb  = MEDIAN INSERT"
          << "  = median insert size; ideal 0.35-0.40" <<  endl;
     PerfStatLoggerM::log( "median_ins_sz", m, "median insert size", true );

     // Print insert behavior stat.

     if ( r2_pct_proper > 0 ) {
          out << "- " << rat(r2_pct_proper,1) << " %   = PROPER PAIRS"
               << "   = fraction of proper read pairs; ideal >=75" <<  endl;
          PerfStatLoggerM::log( 
               "r2_pct_proper", r2_pct_proper, "read two % proper", true );
     }

     // Print molecule length.

     if ( IsRegularFile( DIR + "/patch"+WRITE_SUB+"/a.fhist") ) // temporary!
     {    vec <double> fhist;
          BinaryReader::readFile( DIR + "/patch"+WRITE_SUB+"/a.fhist", &fhist );
          int mol_lwml = 0;
          if ( fhist.size() != 0 ) {
               double l2tot=0.0, ltot=0.0;
               const double BIN_WIDTH=1000.0;
               for (int l = 0; l != fhist.isize(); l++ ) {
                    double bin_mp = (l+0.5)*BIN_WIDTH;
                    l2tot += (fhist[l]*bin_mp*bin_mp);
                    ltot  += (fhist[l]*bin_mp);
               }
               mol_lwml = int( l2tot/ltot );
          }
          out << "- " << rat(mol_lwml,1000) << " kb  = MOLECULE LEN"
               << "   = weighted mean molecule size; ideal 50-100" <<  endl;
          PerfStatLoggerM::log( 
               "lw_mean_mol_len", mol_lwml, 
               "length-weighted mean molecule length", true );    }

     // Print read placement.

     int64_t placed = 0;
     for ( int64_t id = 0; id < (int64_t) dpaths.size( ); id++ )
          if ( dpaths[id].size( ) > 0 ) placed++;
     out << "- " << per(placed,dpaths.size( )) << " %   = PHASED"
          << "         = nonduplicate and phased reads; ideal 45-50" << endl;
     PerfStatLoggerM::log( "placed_fraction", double(placed)/double(dpaths.size( )),
          "nonduplicate and phased reads", true );

     // ============================================================================
     // Print output stats.
     // ============================================================================

     hline( );
     out << "OUTPUT" << endl;
     // out << "- total constituent edges = " << ToStringAddCommas(total_edges)
     //      << endl;
     // PerfStatLoggerM::log( "min_line", MIN_LINE, "minimum line" );
     // out << "- minimum line = " << MIN_LINE << endl;

     // Print scaffold count.

     PerfStatLoggerM::log( "total_lines", nlines/2, "total lines", true);
     out << "- " << rat(nlines/2,1000) << " K   = LONG SCAFFOLDS"
          << " = number of scaffolds >= 10 kb" <<  endl;

     // Print N50 perfect stretch.

     if ( N50_perf >= 0 )
     {    PerfStatLoggerM::log( "N50_perf", N50_perf, "N50 perfect stretch" );
          out << "- " << rat(N50_perf,1000) << " kb  = PERFECT N50"
               << "    = N50 perfect stretch size" <<  endl;    }

     // Print weighted error rate.

     if ( N50_perf >= 0 )
     {    double errwf = int( round( ( errw*1000 ) * 100 ) ) / 100.0;
          PerfStatLoggerM::log( "weighted_error_rate", errwf, "weighted error rate");
          out << "- " << rat(errwf,1) << " /kb = ERROR RATE"
               << "     = weighted error rate" <<  endl;    }

     // Print edge N50.

     out << "- " << rat(N50_edge,1000) << " kb  = EDGE N50" 
          << "       = N50 edge size" <<  endl;
     PerfStatLoggerM::log( "N50_edge", N50_edge, "N50 edge length", true);

     // Print contig N50.

     out << "- " << rat(contig_line_N50,1000) << " kb  = CONTIG N50" 
          << "     = N50 contig size" <<  endl;
     PerfStatLoggerM::log( "N50_contig", contig_line_N50, "N50 contig length", true);

     // Print phase block N50.

     out << "- " << rat(phasetig_N50,1000000) << " Mb  = PHASEBLOCK N50"
          << " = N50 phase block size" << endl;
     PerfStatLoggerM::log( "phasetig_N50", phasetig_N50, "N50 phasetig length", true);

     // Print scaffold N50 and N60.  Treated differently.

     #ifdef CS
     out << "- " << rat(scaffold_N50,1000000) << " Mb  = SCAFFOLD N50" 
          << "   = N50 scaffold size" <<  endl;
     #endif
     PerfStatLoggerM::log( 
          "scaffold_N50", scaffold_N50, "N50 scaffold length", true);
     out << "- " << rat(scaffold_N60,1000000) << " Mb  = SCAFFOLD N60" 
          << "   = N60 scaffold size" <<  endl;
     PerfStatLoggerM::log( 
          "scaffold_N60", scaffold_N60, "N60 scaffold length" );

     // Print approximate genome size.

     out << "- " << rat(est_genome_size,1000000000) << " Gb  = ASSEMBLY SIZE"
          << "  = assembly size (only scaffolds >= 10 kb)" << endl;
     PerfStatLoggerM::log( 
          "assembly_size", est_genome_size, "assembly size", true);

     // Print misassembly stats.

     if ( genome.size( ) > 0 )
     {    PerfStatLoggerM::log( "dis_err_perc",
               nper(total_err_dis_num,total_err_dis_den),
               "fraction of assembly in distant misjoin" );
          out << "- " << per( total_err_dis_num, total_err_dis_den )
               << " %   = dis error"
               << "      = fraction of assembly in distant misjoin" << endl;
          PerfStatLoggerM::log( "ori_err_perc",
               nper(total_err_ori_num, total_err_ori_den ),
               "fraction of assembly having wrong orientation" );
          out << "- " << per( total_err_ori_num, total_err_ori_den )
               << " %   = ori error"
               << "      = fraction of assembly having wrong orientation" << endl;
          PerfStatLoggerM::log( "ord_err_perc",
               nper( total_err_ord_num, total_err_ord_den ),
               "fraction of assembly out of order" );
          out << "- " << per( total_err_ord_num, total_err_ord_den )
               << " %   = ord error"
               << "      = fraction of assembly out of order" << endl;
          double mis = double(total_err_dis_num)/double(total_err_dis_den)
               + double(total_err_ori_num)/double(total_err_ori_den)
               + double(total_err_ord_num)/double(total_err_ord_den);
          PerfStatLoggerM::log( "misasm_perc", mis, "percent misassembled" );
          out << "- " << per( mis*10000, 10000 ) 
               << " %   = MISASSEMBLED" << endl;    }

     // Print gap stats.

     if ( genome.size( ) > 0 )
     {    cap_frac = cap_gap_total / double(genome_size);
          PerfStatLoggerM::log( "air_err_perc", cap_frac,
                    "fraction of assembly in captured gaps");
          out << "- " << per( cap_frac*10000, 10000 )    << " %   = air error"
               << "      = fraction of assembly in captured gaps" << endl;
          double uncap_frac = (genome_size-cov_total) / double(genome_size);
          PerfStatLoggerM::log( "vac_err_perc",
                    uncap_frac, "fraction of assembly in uncaptured gaps");
          out << "- " << per( uncap_frac*10000, 10000 )    << " %   = vac error"
               << "      = fraction of assembly in uncaptured gaps" << endl;
          PerfStatLoggerM::log("total_gap_perc",
                    cap_frac+uncap_frac, "total fraction of assembly in gaps");
          out << "- " << per( (cap_frac+uncap_frac)*10000, 10000 )
               << " %   = GAP" << endl;    }
     hline( );

     // Tidy up.

     // write only CS facing stats
     PerfStatLoggerM::dump_csv( OUTDIR + "/summary_cs.csv" );
     // and all the stats here
     PerfStatLoggerM::dump_csv( OUTDIR + "/all_stats.csv", ",", false );
}
