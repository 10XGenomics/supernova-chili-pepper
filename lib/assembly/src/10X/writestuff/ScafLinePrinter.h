// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// ScafLinePrinter.  Walks and dumps lines of lines from a supergraph.

#ifndef _SCAF_LINE_PRINTER_H
#define _SCAF_LINE_PRINTER_H
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "10X/astats/AssemblyStats.h"
#include "10X/astats/LineLine.h"
#include "10X/DfTools.h"
#include "10X/Gap.h"

#define fwLog(fw, msg) { if ( (fw).IsLog() ) { (fw).Log() << msg; } }

class FastaEdgeWriter {
public:
     FastaEdgeWriter( String filename, String version, int K = 48, Bool abbrev = True,
               Bool gap_absorb = True, String log = "", 
               Bool gap_skip = False, 
               Bool seq_counts = False, unsigned int minsize = 0 ) : 
               _gap_repr_size(100), _K(K), _abbrev(abbrev), 
               _gap_absorb(gap_absorb), 
               _log(log),
               _gap_skip(gap_skip),
               _seq_counts(seq_counts),
               _minsize(minsize),
               _out(filename.c_str() ), _vleft(-1), _vright(-1), _count(0u), _version(version) {
                    if ( _log != "" ) {
                         _log_out.open( log.c_str(), ios::trunc );
                         _islog=True;
                    } else {
                         _log_out.open( "/dev/null" );      // probably some stream sink that we could use instead
                         _islog=False;
                    }

                    _seq.reserve(8000000);
               };

     ~FastaEdgeWriter() { 
          if ( _seq.size() )  {
               cout << "BUG: FastaEdgeWriter: we ended without a Break()" << endl;
          }
     }

     Bool IsLog() { return _islog; }

     void AddGapEdge( int vleft, int vright, int edgeno ) {
          _vright=vright;
          if ( _gap_skip ) {
               fwLog(*this,  "SKIP gap altogether" << endl);
               return;
          }
          if ( _seq.size() > 0 ) {
               PushTail();
          } else {
               ForceAssertEq( _tail.size(), 0u );
               _vleft = vleft;
          }
          if ( !_gap_absorb) {
               this->Break();
               _vleft = vleft;
               fwLog(*this, "YES BREAKING GAP!" << endl);
          } else {
               fwLog(*this, "NOT BREAKING GAP!" << endl);
          }
          _edge_numbers.push_back(edgeno);
          for ( int i = 0; i < _gap_repr_size; i++ )
               _seq.push_back('N');
          if ( !_gap_absorb) {
               this->Break();
          }
     }

     void AddEdge( int vleft, int vright, int edgeno, basevector this_seq ) {
          AddEdgeSeq(vleft, vright, {edgeno}, this_seq);
     }

     void AddEdgeSeq( int vleft, int vright, vec<int> edgenos, basevector this_seq ) {
          _vright=vright;
          _edge_numbers.append(edgenos);
          if ( _seq.size() ) {
               if ( _tail.size() ) {
                    basevector head(this_seq, 0, (_K-1) );
                    if ( head != _tail ) { 
                         Log(" [HEAD!=TAIL] "); 
                         if ( edgenos.size() == 1 ) {
                              Log(" +GAP ");
                              AddGapEdge(vleft, vright, edgenos[0]);
                              return;
                         } else {
                              cout << "OUCH" << endl;
                         }
                    } else { Log(" [head==tail] "); }
               }
          } else {
               _vleft = vleft;
          }

          std::transform( this_seq.begin(), this_seq.end() - (_K-1),
                    std::back_inserter(_seq), BaseToCharMapper() );

          _tail.SetToSubOf( this_seq, this_seq.size() - (_K-1), (_K-1) );
     }

     void PushTail() {
          std::transform(_tail.begin(), _tail.end(),
                    std::back_inserter(_seq), BaseToCharMapper() );

          _tail.clear();
     }


     void Break() {
          ForceAssertGe(_vright, 0);
          fwLog(*this,"      Break for Record " << _count << endl);
          PushTail();

          if ( _seq.size() >= _minsize ) {
               _out << ">";
               if ( _seq_counts) _out << _count << " ";
               _out << "edges=";
               if ( _abbrev ) _out << _edge_numbers.front() << ".." << _edge_numbers.back();
               else _out << printSeq(_edge_numbers);
               _out << " ";
               _out << "left=" << _vleft << " right=" << _vright << " ver=" << _version << endl;
               for ( size_t i = 0; i < _seq.size(); ++i ) {
                    _out << _seq[i];
                    if ( i+1 == _seq.size() || (i%80)==79 )
                         _out << endl;
               }
               _count++;
          } else {
               fwLog(*this, "      -> skipped due to length=" << _seq.size() << endl);
          }

          _seq.clear();
          _seq.reserve(8000000);
          _edge_numbers.clear();
     }

     void BreakIfSeq() {
          if ( _seq.size() > 0 ) {
               Break();
          } else {
               fwLog(*this, "      BreakIfSeq - no sequence" << endl);
          }
     }

     std::ostream& Log( String const& message = "" ) {
          if ( _islog ) return (_log_out << message);
          else return _log_out;
     }

     string Seq() { return _seq; }
     basevector Tail() { return _tail; }

     size_t Count() { return _count; }


private:
     FastaEdgeWriter() = delete;
     FastaEdgeWriter( FastaEdgeWriter const& ) = delete;
     FastaEdgeWriter& operator=( FastaEdgeWriter const& ) = delete;

     int const _gap_repr_size;

     int const _K;
     Bool _abbrev;
     Bool _gap_absorb;
     String _log;
     Bool _islog;
     Bool _gap_skip;
     Bool _seq_counts;
     unsigned int _minsize;

     std::ofstream _out;
     std::ofstream _log_out;

     int _vleft;
     int _vright;
     vec<int> _edge_numbers;
     string _seq;
     basevector _tail;
     size_t _count;

     String _version;
};

class ScafLinePrinter {
public:
     ScafLinePrinter( HyperBasevectorX const& hb, HyperBasevectorX const& hbd, digraphE<vec<int>> const& D,
                        vec<int> const& dinv, LineVec const& dlines,
                        vec<uint16_t> const& dedge_counts, 
                        digraphE<vec<int>> const& D2, HyperBasevectorX const& hbd2 ) : 
                         _hb(hb), _hbd(hbd), _D(D), _dinv(dinv), _dlines(dlines),
                         _dedge_counts( dedge_counts ), 
                         _keepRc(False), _separateGaps(False), _breakBubbles(False),
                         _mashMegaBubbles(False), _seqCounts(False)
     {
          cout << Date() << ": building lines-of-lines" << endl;
          FindLineLines(_D, _dinv, _dlines, _dlines2, &_linv2 );

          D.ToLeft(_to_left);
          D.ToRight(_to_right);
          D2.ToLeft(_to_left2);
          D2.ToRight(_to_right2);

          FindIncorrectSeqInsertions(hbd2, _baddies);
          ReinsertLoopsMap( _D, _dinv, _dcell_map );
     };

     void SetKeepRc( Bool keeprc ) { _keepRc = keeprc; };
     void SetSeparateGaps( Bool sepgap ) { _separateGaps = sepgap; };
     void SetBreakBubbles( Bool breakbubble ) { _breakBubbles = breakbubble; };
     void SetMashMegaBubbles( Bool mashmega ) { _mashMegaBubbles = mashmega; };

     void WalkScaffoldLines( FastaEdgeWriter& fw, size_t choose = 0 ); 

private:
     void FindIncorrectSeqInsertions( HyperBasevectorX const& hbd2, vec<bool>& baddies );
     void ExpandMegabubbleArm( FastaEdgeWriter& ,  vec<int> const& arm );
     int ExpandDLineCell( FastaEdgeWriter&, vec<vec<int>> const& cell ) ;
     void ExpandDLinePath( FastaEdgeWriter&, vec<int> const& path );
     void ExpandDLineGapCell( FastaEdgeWriter& fw, vec<vec<vec<int>>> const& line, int const celli );
     void ExpandDLine( FastaEdgeWriter& fw, int linei );
     void ExpandDGraphGapEdge( FastaEdgeWriter& fw, int edge  );
     void ExpandDGraphCellEdge( FastaEdgeWriter& fw, int edge_id, cell const& this_cell );
     void BustMegabubble( FastaEdgeWriter&, vec<vec<int>> const& cell );

     HyperBasevectorX const& _hb;
     HyperBasevectorX const& _hbd;
     digraphE<vec<int>> const& _D;
     vec<int> const& _dinv;
     LineVec const& _dlines;
     vec<uint16_t> _dedge_counts;

     Bool _keepRc;
     Bool _separateGaps;
     Bool _breakBubbles;
     Bool _mashMegaBubbles;
     Bool _seqCounts;

     LineVec _dlines2;
     vec<int> _linv2;
     vec<int> _to_left;
     vec<int> _to_right;
     vec<int> _to_left2;
     vec<int> _to_right2;
     map<int,int> _dcell_map;
     vec<bool> _baddies;  // bad sequence insertions
};
#endif /* _SCAF_LINE_PRINTER_H */
