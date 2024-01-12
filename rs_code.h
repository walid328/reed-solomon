#ifndef RS_CODE_H
#define RS_CODE_H

#include "polynomial.h"

// An implementation of Reed Solomon codes using Gao's method.

// Encodes the message m = (m_1, m_2, ..., m_k) using the values
// (a_1, a_2, ..., a_n) and stores the associated RS codeworld
// (c_1, c_2, ..., c_n) in c.
void rs_encode(int *c, int *a, int n, int *m, int k);

// Decodes the received code word b and stores the message
// m. If there is too many errors quits and displays and
// error message.
void rs_decode(int *m, int *a, int *b, int n, int k);

#endif