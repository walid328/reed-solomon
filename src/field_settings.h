#ifndef FIELD_SETTINGS_H
#define FIELD_SETTINGS_H

#include "finite_field.h"
#include "zp_array.h"

// p is the size of the field.
// p < 2^31 to avoid overflow.
// <==> p <= 2147483647
// p = q*n + 1 where n is a power of 2 and q is odd.
// omega is a n^th root of unity.
// omegas is an zp_array containing the powers of omega.

extern int p;
extern int q;
extern int n;
extern zp_t omega;
extern zp_array omegas;
extern zp_array inverses;

void field_settings_set(int p_);

void field_settings_set(int p_);

void field_settings_free(void);

#endif
