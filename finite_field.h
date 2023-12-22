#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H
#define p 307

// Basic operations in Z/pZ.

// Returns the sum of two elements of Z/pZ.
int add_zp(int n, int m);

// Returns the difference of two elements of Z/pZ.
int sub_zp(int n, int m);

// Returns the opposite of an element of Z/pZ.
int opp_zp(int n);

// Returns the product of two elements of Z/pZ.
int mul_zp(int n, int m);

// Return the inverse of n in Z/pZ.
int inv_zp(int n);

// Returns a random element in Z/pZ.
int rand_zp();

#endif