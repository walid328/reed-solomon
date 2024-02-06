#ifndef FIELD_SETTINGS_H
#define FIELD_SETTINGS_H

#include "array.h"

// p is the size of the field.
// (p-1) + (p-1) <= 2^31 - 1 to avoid overflow.
// <==> p <= 2147483648
// p = q*n + 1 where n is a power of 2 and q is odd.
// omega is a n^th root of unity.
// omegas is an array containing the powers of omega.

extern int p;
extern int q;
extern int n;
extern int omega;
extern array omegas;
extern array inverses;

void field_settings_set(int p_);

void field_settings_free(void);

#endif
