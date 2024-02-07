#ifndef RS_CODE_H
#define RS_CODE_H

#include "array.h"

// An implementation of Reed Solomon codes using Gao's method
// for decoding.

/******************************************************/

/* Encoding functions */

// Encode the message using RS(block_length, message_length) and return the
// corresponding codeword. message should have a size
// of message_length at most. points should have block_length distinct elements.
array rs_encode(unsigned int block_length, unsigned int message_length, array points, array message);

// Same as above but points can't be choosen and should be powers of a
// primitive block_length^th root. block_length should be a power of 2 smaller than n.
array rs_encode_2(unsigned int block_length, unsigned int message_length, array message);

// Same as above but use fast operations.
array rs_fast_encode(unsigned int block_length, unsigned int message_length, array message);

/******************************************************/

/* Decoding function */

// Decode the received message and return the corresponding
// message. received should have a size of block_length at most. points
// should have block_length distinct elements.
array rs_decode(unsigned int block_length, unsigned int message_length, array points, array received);

// Same as above but points can't be choosen and should be powers of a
// primitive block_length^th root. block_length should be a power of 2 smaller than n.
array rs_decode_2(unsigned int block_length, unsigned int message_length, array received);

// Same as above but use fast operations.
array rs_fast_decode(unsigned int block_length, unsigned int message_length, array received);

#endif