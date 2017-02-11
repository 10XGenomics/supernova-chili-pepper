//
//  Author: Preyas Shah - Apr 08, 2016 - 10x
//

#include "10X/paths/ReadPathVecX.h"

#define READPATHVECX_CC_
#ifdef READPATHVECX_CC_

// A description of a graph traversal by some sequence (a read, let's say).
// The logical representation of this class in memory is in a {numEdges,CharArray} where
// the data is packed serially as:
// {
// uchar numEdges
// int16_t offset
// uint32_t Edge1
// branch2-branch(N-1) as 2bit representation OR other coding scheme
// }
// This is a major revision of an existing sturcture "ReadPath" (under /paths/long/) that allows efficient storage of data on disk and in memory at a low cost to random access and unpack.

// TODO
// make thread-safe?
// throw away useless functions later
// provide binary r/w functionality like the feudal system
// Redo for HyperBasevector?
// the implementation of jump ptr will change if we allow more than 4 branch ids

/* I T E R A T O R S */

/**
 * @brief return begin itr to internal storage
 *
 * @return 
 */
uchararr_iterator ReadPathVecX::begin() { return ZippedData.begin(); }

/**
 * @brief return end itr to internal storage
 *
 * @return 
 */
uchararr_iterator ReadPathVecX::end() { return ZippedData.end(); }

/**
 * @brief return begin const itr
 *
 * @return 
 */
const_uchararr_iterator ReadPathVecX::begin() const { return ZippedData.begin(); }

/**
 * @brief return end const itr
 *
 * @return 
 */
const_uchararr_iterator ReadPathVecX::end() const { return ZippedData.end(); }

/**
 * @brief return cbegin const itr
 *
 * @return 
 */
const_uchararr_iterator ReadPathVecX::cbegin() const { return ZippedData.cbegin(); }

/**
 * @brief return cend const itr
 *
 * @return 
 */
const_uchararr_iterator ReadPathVecX::cend() const { return ZippedData.cend(); }

/**
 * @brief iterator to a random read inside the compressed storage
 *
 * @param read_id: const int64_t
 *
 * @return 
 */
const_uchararr_iterator ReadPathVecX::accessItr(const int64_t read_id) const{
    ForceAssertLt(read_id,next_start_rid);
    ForceAssertGt(read_id+1, start_rid);
    int64_t floor_idx = (read_id-start_rid)/skip;

    const_uchararr_iterator itr = ZippedData.begin() + ZipIndex[floor_idx];
    for(int64_t i = start_rid+floor_idx*skip; i<read_id; i++){
        jumpPtr(itr);
    }
    return itr;
}

/**
 * @brief jumps an input ptr to next read starting byte, use sparingly
 *
 * @param itr
 */
void ReadPathVecX::jumpPtr(const_uchararr_iterator& itr) const{
    unsigned char numEdges = *itr;
    // this would change if branch id not standard
    itr = itr + ( (static_cast<int>(numEdges)>0)?1+2+4+((numEdges-1)+3)/4:1 );
}

/**
 * @brief index to a random read inside the compressed storage
 *
 * @param read_id: const int64_t
 *
 * @return 
 */
int64_t ReadPathVecX::accessIdx(const int64_t read_id) const{
    ForceAssertLt(read_id,next_start_rid);
    ForceAssertGt(read_id+1, start_rid);
    int64_t floor_idx = (read_id-start_rid)/skip;

    int64_t idx = ZipIndex[floor_idx];
    for(int64_t i = start_rid+floor_idx*skip; i<read_id; i++){
        jumpIdx(idx);
    }
    return idx;
}

/**
 * @brief provides next readid index in storage, use sparingly
 *
 * @param idx: int64_t&
 */
void ReadPathVecX::jumpIdx(int64_t& idx) const{
    unsigned char numEdges = ZippedData[idx];
    // this would change if branch id not standard
    idx += ( (static_cast<int>(numEdges)>0)?1+2+4+((numEdges-1)+3)/4:1 );
}


/* C O N S T R U C T O R S */

/**
 * @brief default constructor
 */
ReadPathVecX::ReadPathVecX(){}

/**
 * @brief reserve appropriate space using heuristics. expect ~ 10byte meanlength we expect at least a few thousand reads per RPY
 *
 * @param numReads: const int64_t
 */
ReadPathVecX::ReadPathVecX(const int64_t numReads){
    reserve(numReads);
}

/**
 * @brief reserve appropriate space using heuristics, from sid to sid+numReads
 *
 * @param sid: const int64_t
 * @param numReads: const int64_t
 */
ReadPathVecX::ReadPathVecX(const int64_t sid, const int64_t numReads){
    reserve(sid,numReads);
}

/**
 * @brief appends the zipped sequence to internal storage
 *
 * @param rp: const ReadPath&
 * @param hb: const HyperBasevectorX&
 */
ReadPathVecX::ReadPathVecX(const ReadPath& rp, const HyperBasevectorX& hb){
    append(rp,hb);
}

/**
 * @brief same as above, accepts vec<int> and offset as argument explicitly
 *
 * @param edge_list: const vec<int>&
 * @param offset: const int
 * @param hb: const HyperBasevectorX&
 */
ReadPathVecX::ReadPathVecX(const vec<int>& edge_list, const int offset, const HyperBasevectorX& hb){
    append(edge_list,offset,hb);
}

/**
 * @brief constructor takes and appends a ReadPathX
 *
 * @param rpx: const ReadPathX
 */
ReadPathVecX::ReadPathVecX(const ReadPathX rpx){
    append(rpx);
}

/**
 * @brief constructor takes a MasterVec of ReadPaths and appends it
 *
 * @param paths: const ReadPathVec&
 * @param hb: const HyperBasevectorX&
 */
ReadPathVecX::ReadPathVecX(const ReadPathVec& paths, const HyperBasevectorX& hb){
    append(paths, hb);
}

/**
 * @brief constructor takes a MasterVec of ReadPaths and appends from sid to sid+numReads
 *
 * @param paths: const ReadPathVec&
 * @param hb: const HyperBasevectorX&
 * @param sid: const int64_t
 * @param numReads: const int64_t
 */
ReadPathVecX::ReadPathVecX(const ReadPathVec& paths, const HyperBasevectorX& hb, const int64_t sid, const int64_t numReads){
    append(paths, hb, sid, numReads);
}

/**
 * @brief constructor takes a virtual MasterVec of ReadPaths and appends
 *
 * @param paths: const VirtualMasterVec<ReadPath>&
 * @param hb: const HyperBasevectorX&
 */
ReadPathVecX::ReadPathVecX(const VReadPathVec& paths, const HyperBasevectorX& hb){
    append(paths, hb);
}

/**
 * @brief constructor takes a virtual MasterVec of ReadPaths and appends from sid to sid+numReads
 *
 * @param paths: const VirtualMasterVec<ReadPath>&
 * @param hb: const HyperBasevectorX&
 * @param sid: const int64_t
 * @param numReads: const int64_t
 */
ReadPathVecX::ReadPathVecX(const VReadPathVec& paths, const HyperBasevectorX& hb, const int64_t sid, const int64_t numReads){
    append(paths, hb, sid, numReads);
}


/* P S E U D O - C O N S T R U C T O R S */

/**
 * @brief helper, reserve space worth numReads and update internal state params
 *
 * @param numReads: 
 */
void ReadPathVecX::reserve(const int64_t numReads){
    ForceAssertGt(numReads,0);
    start_rid = 0;
    next_start_rid = 0;
    ZipIndex.reserve(numReads/skip+3);
    ZipIndex.push_back(0);
    int64_t reserve_bytes = static_cast<int64_t>(numReads*10);
    ZippedData.reserve(reserve_bytes);
}

/**
 * @brief helper, reserve space given sid and numReads
 *
 * @param sid
 * @param numReads
 */
void ReadPathVecX::reserve(const int64_t sid, const int64_t numReads){
    ForceAssertGt(numReads,0);
    start_rid = sid;
    next_start_rid = sid;
    ZipIndex.reserve(numReads/skip+3);
    ZipIndex.push_back(0);
    int64_t reserve_bytes = numReads*10;
    ZippedData.reserve(reserve_bytes);
}

/**
 * @brief appends ReadPath; ideally should reserve mem calling 2nd constructor prior
 *
 * @param rp: const ReadPath&
 * @param hb: const HyperBasevectorX&
 */
void ReadPathVecX::append(const ReadPath& rp, const HyperBasevectorX& hb){
    LLzip(ZippedData,rp,hb);
    updateZipIndex();
}

/**
 * @brief same as above, accepts vec<int> and offset as argument explicitly
 *
 * @param edge_list: const vec<int>&
 * @param offset: const int
 * @param hb: const HyperBasevectorX&
 */
void ReadPathVecX::append(const vec<int>& edge_list, const int offset, const HyperBasevectorX& hb){
    LLzip(ZippedData,edge_list,offset,hb);
    updateZipIndex();
}

/**
 * @brief appends raw compressed seq
 *
 * @param rpx: const ReadPathX&
 */
void ReadPathVecX::append(const ReadPathX& rpx){
    ZippedData.insert(ZippedData.end(), rpx.begin(), rpx.end());
    updateZipIndex();
}

/**
 * @brief appends raw RPY data; do not pass self argument!
 *
 * @param rpvx: const ReadPathVecX&
 */
void ReadPathVecX::append(const ReadPathVecX& rpvx){
    if(&rpvx==this){
        cout<<"Do not pass self as argument!"<<endl;
        exit(1);
    }
    ZippedData.reserve(storageSize()+rpvx.storageSize()*10);//heuristic
    ZippedData.insert(ZippedData.end(),rpvx.begin(),rpvx.end());
    recompileZipIndex();
}

/**
 * @brief appends MasterVec<ReadPath>
 *
 * @param paths: const ReadPathVec&
 * @param hb: const HyperBasevectorX&
 */
void ReadPathVecX::append(const ReadPathVec& paths, const HyperBasevectorX& hb){
    reserve(paths.size());
    for(uint64_t id = 0; id< paths.size(); id++){
        append(paths[id],hb);
    }
}

/**
 * @brief appends MasterVec<ReadPath> from sid to sid+numReads
 *
 * @param paths: const ReadPathVec&
 * @param hb: const HyperBasevectorX&
 * @param sid: const int64_t
 * @param numReads: const int64_t
 */
void ReadPathVecX::append(const ReadPathVec& paths, const HyperBasevectorX& hb, const int64_t sid, const int64_t numReads){
    reserve(sid,numReads);
    for(int64_t id = sid; id< sid+numReads; id++){
        append(paths[id],hb);
    }
}

/**
 * @brief appends VirtualMasterVec<ReadPaths>
 *
 * @param paths: const VirtualMasterVec<ReadPath>&
 * @param hb: const HyperBasevectorX&
 */
void ReadPathVecX::append(const VReadPathVec& paths, const HyperBasevectorX& hb){
    reserve(paths.size());
    for(uint64_t id = 0; id< paths.size(); id++){
        append(paths[id],hb);
    }
}

/**
 * @brief appends VirtualMV<ReadPaths> from sid to sid+numReads
 *
 * @param paths: const VirtualMasterVec<ReadPath>&
 * @param hb: const HyperBasevectorX&
 * @param sid: const int64_t
 * @param numReads: const int64_t
 */
void ReadPathVecX::append(const VReadPathVec& paths, const HyperBasevectorX& hb, const int64_t sid, const int64_t numReads){
    reserve(sid,numReads);
    for(int64_t id = sid; id< sid+numReads; id++){
        append(paths[id],hb);
    }
}

/**
 * @brief unpacks a ReadPathVecX into a ReadPath
 *
 * @param rp: ReadPath&
 * @param hb: const HyperBasevectorX&
 * @param read_id: const int64_t
 */
void ReadPathVecX::unzip(ReadPath& rp, const HyperBasevectorX& hb, const int64_t read_id) const{
    LLunzip(ZippedData,rp,hb,accessIdx(read_id));
}

/**
 * @brief unpacks a ReadPathVecX into an {edge list, offset} for given read id
 *
 * @param edge_list: vec<int>&
 * @param offset: int&
 * @param hb: const HyperBasevectorX&
 * @param read_id: const int64_t
 */
void ReadPathVecX::unzip(vec<int>& edge_list, int& offset, const HyperBasevectorX& hb, const int64_t read_id) const {
    LLunzip(ZippedData,edge_list,offset,hb,accessIdx(read_id));
}

/**
 * @brief unzips internal data into a MasterVec of ReadPaths
 *
 * @param paths: ReadPathVec&
 * @param hb: const HyperBasevectorX&
 */
void ReadPathVecX::unzip(ReadPathVec& paths, const HyperBasevectorX& hb) const{
    // clear existing paths
    paths.clear();
    paths.reserve(size());
    for(int64_t id = start_rid; id<next_start_rid; id++){
        ReadPath rp;
        unzip(rp,hb,id);
        paths.push_back(rp);
    }
}

/**
 * @brief unzips internal data and appends to existing MasterVec of ReadPaths, quite slow
 *
 * @param paths: ReadPathVec&
 * @param hb: const HyperBasevectorX&
 */
void ReadPathVecX::unzipAppend(ReadPathVec& paths, const HyperBasevectorX& hb) const{
    paths.reserve(paths.size()+size());
    for(int64_t id = start_rid; id<next_start_rid; id++){
        ReadPath rp;
        unzip(rp,hb,id);
        paths.push_back(rp);
    }
}


/* C O P Y   &   A S S I G N M E N T */
/**
 * @brief deep copy constr
 *
 * @param source: const ReadPathVecX&
 */
ReadPathVecX::ReadPathVecX(const ReadPathVecX& source){
    skip = source.skip;
    start_rid = source.start_rid;
    next_start_rid = source.next_start_rid;
    ZippedData = source.ZippedData;
    ZipIndex = source.ZipIndex;
}

/**
 * @brief operator=
 *
 * @param source: const ReadPathVecX&
 *
 * @return 
 */
ReadPathVecX& ReadPathVecX::operator=(const ReadPathVecX& source){
    if (this == &source)
        return *this;
    skip = source.skip;
    start_rid = source.start_rid;
    next_start_rid = source.next_start_rid;
    ZippedData = source.ZippedData;
    ZipIndex = source.ZipIndex;
    return *this;
}
/* U P D A T E   &   M A I N T A I N   I N D E X */

/**
 * @brief updates the ZipIndex every skip read_ids, updates the next_start_rid 
 */
void ReadPathVecX::updateZipIndex(){
    // this will be fast, assuming reserved, or at least initialized
    // store every skip^th read_id, starting from start_rid
    next_start_rid++;
    if ((next_start_rid-start_rid)%skip == 0){
        ZipIndex.push_back(storageSize());
    }
}

/**
 * @brief dangerous private functions; depends on accessIdx() not accessing ZipIndex out of bounds
 */
void ReadPathVecX::recompileZipIndex(){
    int64_t rid = skip*ZipIndex.size()+start_rid+1;
    ZipIndex.reserve(storageSize()/meanPathLength()/skip+5);
    int64_t accidx = accessIdx(rid);
    while(accidx<storageSize()){
        rid++;
        accidx = accessIdx(rid);
        if((rid-start_rid)%skip==0)
            ZipIndex.push_back(accidx);
    }
    next_start_rid = ZipIndex.size()*skip + start_rid;
}

/* D E S T R U C T O R */

/**
 * @brief destructor
 */
ReadPathVecX::~ReadPathVecX(){
    ZippedData.clear();
    ZipIndex.clear();
}

void ReadPathVecX::Destroy(){
    start_rid = 0;
    next_start_rid = 0;
    skip = 0;
    ZippedData.clear();
    ZipIndex.clear();
}


/* P S E U D O   M O D I F I E R S   &   A C C E S S O R S  F U N C T I O N S */

/**
 * @brief returns total interal storage size in bytes
 *
 * @return 
 */
int64_t ReadPathVecX::storageSize() const{
    return ZippedData.size();
}

/**
 * @brief return number of reads
 *
 * @return 
 */
int64_t ReadPathVecX::size() const{
    return (next_start_rid-start_rid);
}

/**
 * @brief returns number of reads
 *
 * @return 
 */
int64_t ReadPathVecX::numReadPaths() const{
    return (next_start_rid-start_rid);
}

/**
 * @brief returns mean storage space per read
 *
 * @return 
 */
double ReadPathVecX::meanPathLength() const{
    return static_cast<double>(storageSize())/(next_start_rid-start_rid);
}

/**
 * @brief extracts a rpx subvector corresponding to a read id
 *
 * @param rpx: ReadPathX&
 * @param read_id: const int64_t
 */
void ReadPathVecX::extractRPX(ReadPathX& rpx, const int64_t read_id){
    rpx.assign(accessItr(read_id),accessItr(read_id+1));
}

/**
 * @brief accessors for numEdges for given read id
 *
 * @param read_id: const int64_t
 *
 * @return 
 */
unsigned char ReadPathVecX::getNumEdges(const int64_t read_id) const {
    return LLgetNumEdges(ZippedData,accessIdx(read_id));
}

/**
 * @brief accessors for offsets at given read_id
 *
 * @param read_id: const int64_t
 *
 * @return 
 */
int ReadPathVecX::getOffset(const int64_t read_id) const { 
    return LLgetOffset(ZippedData,accessIdx(read_id));
}

/**
 * @brief sets offset for given read id
 *
 * @param read_id: const int64_t
 * @param offset: int
 */
void ReadPathVecX::setOffset(const int64_t read_id, int offset) { 
    // TODO 
    FatalErr("setOffset is not implemented yet");
}

/**
 * @brief add offset at given read_id
 *
 * @param read_id: const int64_t
 * @param add: int
 */
void ReadPathVecX::addOffset(const int64_t read_id, int add ) {
    // TODO
    FatalErr("addOffset is not implemented yet");
}

/**
 * @brief gets first skip for given read id
 *
 * @param read_id: const int64_t
 *
 * @return 
 */
unsigned ReadPathVecX::getFirstSkip(const int64_t read_id) const { 
    int mOffset = LLgetOffset(ZippedData,accessIdx(read_id));
    return (mOffset < 0 ? 0u : static_cast<unsigned>(mOffset));
}

/**
 * @brief sets first skip for given read id
 *
 * @param read_id: const int64_t
 * @param firstSkip: unsigned
 */
void ReadPathVecX::setFirstSkip(const int64_t read_id, unsigned firstSkip ) { 
    // TODO
    FatalErr("setFirstSkip is not implemented yet");
}


/* S T D   O/P */

/**
 * @brief outputs unformatted bytes in internal storage; You don't want to use this- outputs too much stuff
 *
 * @param os: ostream&
 * @param rpvx: const ReadPathVecX&
 *
 * @return 
 */
ostream& operator<<( ostream& os, const ReadPathVecX& rpvx)
{
    // std::copy( rpvx.ZippedData.begin(), rpvx.ZippedData.end(), std::ostream_iterator<unsigned int>(os, " ") );
    cout<<"Zipped Data: [";
    for(int64_t i = 0; i<rpvx.storageSize(); i++)
        cout<<static_cast<int>(rpvx.ZippedData[i])<<",";
    cout<<"]"<<endl<<"Zip Index: [";
    for(uint64_t i = 0; i<rpvx.ZipIndex.size(); i++)
        cout<<static_cast<int>(rpvx.ZipIndex[i])<<",";
    cout<<"]"<<endl<<"start_rid: "<<rpvx.start_rid<<", next_start_rid: "<<rpvx.next_start_rid<<", skip: "<<rpvx.skip;
    return os;
}

/**
 * @brief prints a well formatted readpath
 *
 * @param hb: HyperBasevectorX&
 * @param messg: const std::string, default ""
 * @param index: const int64_t, default 0
 */
void ReadPathVecX::printReadPath(const HyperBasevectorX& hb, const std::string messg, const int64_t index) const{
    LLprintMe(ZippedData,hb,messg,index);
}


/* R E A D  -  W R I T E */

/**
 * @brief write binary: SKIP|STARTRID|NEXTSTARTRID|ZIPINDEX|ZIPPEDDATA
 *
 * @param pathname: std::string
 */
void ReadPathVecX::readBinary(std::string pathname){


    std::ifstream INFILE(pathname, std::ios::in | std::ifstream::binary);
    if(!INFILE.is_open()){
        cerr<<"Error opening paths file at "<<pathname<<endl;
        exit(1);
    }

    // read skip, start_rid, next_start_rid
    INFILE.read(reinterpret_cast<char *>(&skip), sizeof(skip));
    INFILE.read(reinterpret_cast<char *>(&start_rid), sizeof(start_rid));
    INFILE.read(reinterpret_cast<char *>(&next_start_rid), sizeof(next_start_rid));
    int64_t ziSize;
    INFILE.read(reinterpret_cast<char *>(&ziSize), sizeof(ziSize));
    int64_t zdSize;
    INFILE.read(reinterpret_cast<char *>(&zdSize), sizeof(zdSize));

    // read ZipIndex
    ZipIndex.clear();
    ZipIndex.reserve(ziSize+3);;
    int64_t zi;
    for(int64_t i = 0; i<ziSize; i++){
        INFILE.read(reinterpret_cast<char *>(&zi), sizeof(zi));
        ZipIndex.push_back(zi);
    }

    // read ZippedData
    ZippedData.clear();
    ZippedData.reserve(zdSize+4);
    unsigned char zd;
    for(int64_t i = 0; i<zdSize; i++){
        INFILE.read(reinterpret_cast<char *>(&zd), sizeof(zd));
        ZippedData.push_back(zd);
    }
    INFILE.close();

}

/**
 * @brief read binary
 *
 * @param pathname: std::string
 */
void ReadPathVecX::writeBinary(std::string pathname) const{
   
    std::ofstream OUTFILE(pathname, std::ios::out | std::ofstream::binary);

    // write skip, start_rid, next_start_rid
    OUTFILE.write(reinterpret_cast<const char *>(&skip), sizeof(skip));
    OUTFILE.write(reinterpret_cast<const char *>(&start_rid), sizeof(start_rid));
    OUTFILE.write(reinterpret_cast<const char *>(&next_start_rid), sizeof(next_start_rid));
    int64_t ziSize = ZipIndex.size();
    OUTFILE.write(reinterpret_cast<const char *>(&ziSize), sizeof(ziSize));
    int64_t zdSize = ZippedData.size();
    OUTFILE.write(reinterpret_cast<const char *>(&zdSize), sizeof(zdSize));

    // write ZipIndex
    OUTFILE.write(reinterpret_cast<const char *>(ZipIndex.data()), sizeof(int64_t)*ziSize);

    // write ZippedData
    OUTFILE.write(reinterpret_cast<const char *>(ZippedData.data()), sizeof(unsigned char)*zdSize);
    OUTFILE.close();

}

#endif
/* READPATHVECX_CC_ */
