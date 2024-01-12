#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H
#define p 193

// Basic operations in Z/pZ.

// Returns the representant of n in Z/pZ.
int zp_mod(int n);

// Returns the sum of two elements of Z/pZ.
int zp_add(int n, int m);

// Returns the difference of two elements of Z/pZ.
int zp_sub(int n, int m);

// Returns the opposite of an element of Z/pZ.
int zp_opp(int n);

// Returns the product of two elements of Z/pZ.
int zp_mul(int n, int m);

// Return the inverse of n in Z/pZ.
int zp_inv(int n);

// Returns a random element in Z/pZ.
int zp_rand();

// Returns base^exp mod p using modular exponentiation.
int zp_exp(int base, int exp);

// Returns a d^th primitive root of unity in Z/pZ.
// We should have p = dq + 1 where d is a power of 2,
// and q is odd.
int zp_prim_root(int q, int d);

// Returns the minimal d^th primitive root of unity in Z/pZ.
// We should have p = dq + 1 where d is a power of 2,
// and q is odd.
int zp_prim_root_min(int q, int d);

#endif