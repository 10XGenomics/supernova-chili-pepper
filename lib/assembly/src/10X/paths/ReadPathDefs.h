/*
 * ReadPathDefs.h
 *
 *  Created on: Apr 06, 2016
 *      Author: preyas shah
 */

#ifndef READPATHDEFS_H_
#define READPATHDEFS_H_

#include "Vec.h"
#include "paths/long/ReadPath.h"
#include "paths/HyperBasevector.h"

enum CODING_SCHEME{FIXED_WIDTH, HUFFMAN};

typedef vec<unsigned char> uCharArray;

typedef uCharArray::iterator uchararr_iterator;
typedef uCharArray::const_iterator const_uchararr_iterator;

typedef VirtualMasterVec<ReadPath> VReadPathVec;

#endif /* READPATHDEFS_H_ */
