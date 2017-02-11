// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_DF_TOOLS_H
#define TENX_DF_TOOLS_H

#if defined(SHOWMEM)
#define MEM(X) { cout << #X << ": mem = " << MemUsageGBString( ) << ", peak = " << PeakMemUsageGBString( ) << endl; }
#else
#define MEM(X)
#endif

#include "CoreTools.h"
#include "Intvector.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

enum class ReadDataType : uint8_t {
     PCR_FREE,            /* PCR-free conventional data */
     PCR,                 /* PCR-amplified conventional data */
     UNBAR_10X,           /* 10X nominally-barcoded data WITHOUT a passing barcode */
     BAR_10X              /* 10X barcoded data WITH a passing barcode */
};

struct DataSet {
     ReadDataType dt;
     int64_t start;

     friend ostream& operator<<( ostream& out, DataSet const& val ) {
          if ( val.dt == ReadDataType::PCR_FREE ) out << "PCR_FREE";
          else if ( val.dt == ReadDataType::PCR ) out << "PCR";
          else if ( val.dt == ReadDataType::UNBAR_10X ) out << "UNBAR_10X";
          else if ( val.dt == ReadDataType::BAR_10X ) out << "BAR_10X";
          else FatalErr("bad value for ReadDataType");
          out << " starts at " << val.start;
          return out;
     }
};

TRIVIALLY_SERIALIZABLE(DataSet);

class PerfStatLoggerM
{
public:
     static void init( bool silent = true ) {
          delete gInst.mpLines;
          gInst.mpLines = new vec<_line>;
          gInst.mSilent = silent;
     }

     template <class T>
     static void log( String const& statName, T const& val, String const& gloss, bool cs = false )
     {
          if ( gInst.mpLines ) {
               _line tmp{ cs, statName, ToString(val), gloss };
               gInst.mpLines->emplace_back( tmp );

               if ( !gInst.mSilent ) {
                    cout << statName << "=" << val << "\t" << gloss << endl;
               }
          }
     }
     // this function dumps out all statistics into a txt file
     // independent of the value of the flag cs
     static void dump_text( String const& filename ) {
          if ( gInst.mpLines ) {
               ofstream out( filename );
               if ( out ) {
                    for ( auto const& line : *gInst.mpLines ) {
                         out << line.name << "\t" << line.sval << "\t" << line.gloss << endl;
                    }
               } else FatalErr("failed writing output to "+filename );
          }
     }
     // this function is used to define a user-facing csv
     // cs_only = true means only CS facing stats are written
     static void dump_csv( String const& filename, String const sep=",", bool cs_only = true ) {
          if ( gInst.mpLines ) {
               ofstream out( filename );
               if ( out ) {
                    const int num_stats = gInst.mpLines->size();
                    bool first_entry=true;
                    for ( auto const& line : *gInst.mpLines ) {
                         if ( cs_only && !line.cs )
                              continue;
                         if ( !first_entry ) {
                              out << sep;
                         }
                         first_entry=false;
                         out << line.name;
                    }
                    out <<endl;
                    first_entry=true;
                    for ( auto const& line : *gInst.mpLines ) {
                         if ( cs_only && !line.cs )
                              continue;
                         if ( !first_entry ) {
                              out << sep;
                         }
                         first_entry=false;
                         out << line.sval;
                    }
                    out <<endl;
               } else FatalErr("failed writing output to "+filename );
          }
     }
     // TODO: we'll make this polymorphic so it correctly dumps out types
     // I just want the API defined right now.
     static void dump_json( String const& filename ) {
          if ( gInst.mpLines ) {
               ofstream out( filename );
               if ( out ) {
                    out << "{" << endl;
                    for ( auto itr = gInst.mpLines->begin(); itr != gInst.mpLines->end();  ) {
                         auto const& line = *itr;
                         out << '\t' << '"' << line.name << '"' << ": { \"value\": \"" <<
                              line.sval << "\", \"descr\": \"" << line.gloss << "\" }";
                         if ( ++itr != gInst.mpLines->end() ) out << "," << endl;
                    }
                    out << endl << "}" << endl;
               } else FatalErr("failed writing output to " + filename );
          }
     }



private:
     struct _line {
          bool   cs;          // should this stat be displayed in CS
          String name;        // stat name
          String sval;        // stat value
          String gloss;       // description of stat
     };
     vec<_line> *mpLines;
     static PerfStatLoggerM gInst;
     bool mSilent;
};

template<int K> void MapClosures( const HyperBasevectorX& hb,
     const vec<basevector>& closures, ReadPathVec& clop );

vec<int> GetBarcodes( const int e, const vec<int>& inv,
     const VecULongVec& paths_index, const vec<int>& bc );

void LoadData( const String& work_dir, const String& R, const vec<String>& lr,
     const vec<double>& LR_SELECT_FRAC, vecbasevector& bases,
     ObjectManager<VecPQVec>& quals_om, vec<int64_t>& bci,
     vec<String>& subsam_names, vec<int64_t>& subsam_starts, vec<DataSet>& datasets);

void GetQualStats( 
     const VecPQVec& quals, vec<vec<vec<int64_t>>> & hist, int & max_read_length);

void FragDist( const HyperBasevectorX& hb, const vec<int>& inv,
     const ReadPathVec& paths, vec<int64_t>& count );

void ReadTwoPctProper( const HyperBasevectorX& hb, const vec<int>& inv,
     VirtualMasterVec<ReadPath> const& paths, double & r2_pct_proper );

typedef triple<int,int,int> trip;
extern template class SmallVec<trip,MempoolAllocator<trip> >;
extern template class OuterVec< SerfVec<trip> >;
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
template class SmallVec<trip,MempoolAllocator<trip> >;
template class OuterVec< SerfVec<trip> >;

void SanityCheckBarcodeCounts( const vec<int64_t>& bci );

void MakeDots( int& done, int& ndots, const int total );

void SuperToSeqGraph( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     HyperBasevectorX& hbd );

#endif
