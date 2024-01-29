#ifndef FIELD_SETTINGS_H
#define FIELD_SETTINGS_H

#include "array.h"

// p is the size of the field.
// p = q*n + 1 where n is a power of 2 and q is odd.
// omega is a n^th root of unity.
// omegas is an array containing the powers of omega.

extern int p;
extern int q;
extern int n;
extern int omega;
extern array omegas;

void field_settings_update(void);

void field_settings_free(void);

#endif
