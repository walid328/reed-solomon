#ifndef RS_CODE_H
#define RS_CODE_H

#include "zp_array.h"
#include "polynomial.h"

// An implementation of Reed Solomon codes using Gao's method
// for decoding.

/******************************************************/

/* Functions to compute g_0 */

// Computes the polynomial g_0 = \prod_{i = 1}^n (x - a_i)
// where {a_0, ..., a_n} = points, n = block_length.
void rs_g_0(poly g_0, zp_array points, int block_length);

// Computes the polynomial g_0 = \prod_{i = 1}^n (x - a_i)
// where n = block_lengt, {a_0, ..., a_n} are powers of a n^th
// primitive root of unity.
void rs_g_0_fourier(poly g_0, int block_length);

/******************************************************/

/* Encoding functions */

// Encode the message using RS(block_length, message_length) and return the
// corresponding codeword. message should have a size
// of message_length at most. points should have block_length distinct elements.
zp_array rs_encode(int block_length, int message_length, zp_array points, zp_array message);

// Same as above but points can't be choosen and should be powers of a
// primitive block_length^th root. block_length should be a power of 2 smaller than n.
zp_array rs_encode_2(int block_length, int message_length, zp_array message);

// Same as above but use fast operations.
zp_array rs_fast_encode(int block_length, int message_length, zp_array message);

/******************************************************/

/* Decoding function */

// Decode the received message and return the corresponding
// message. received should have a size of block_length at most. points
// should have block_length distinct elements.
zp_array rs_decode(poly g_0, int block_length, int message_length, zp_array points, zp_array received);

// Same as above but points can't be choosen and should be powers of a
// primitive block_length^th root. block_length should be a power of 2 smaller than n.
zp_array rs_decode_2(poly g_0, int block_length, int message_length, zp_array received);

// Same as above but use fast operations.
zp_array rs_fast_decode(poly g_0, int block_length, int message_length, zp_array received);

#endif