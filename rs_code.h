#ifndef RS_CODE_H
#define RS_CODE_H

#include "array.h"

// An implementation of Reed Solomon codes using Gao's method
// for decoding.

// Encode the message using RS(n, k) and return the
// corresponding codeword. message should have a size
// of k at most. points should have n distinct elements.
array rs_encode(int n, int k, array points, array message);

// Decode the received message and return the corresponding
// message. received should have a size of n at most. points
// should have n distinct elements.
array rs_decode(int n, int k, array points, array received);

#endif