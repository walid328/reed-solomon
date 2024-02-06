
#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

// Basic operations in Z/pZ.

// Compute a mod p.
int zp_mod(int a);

// Compute a + b mod p.
int zp_add(int a, int b);

// Compute a - b mod p.
int zp_sub(int a, int b);

// Compute -a mod p.
int zp_opp(int a);

// Compute a*b mod p.
int zp_mul(int a, int b);

// Compute a^(-1) mod p via xgcd.
int zp_inv_xgcd(int a);

// Return a^(-1) mod p via precomputed table.
int zp_inv(int a);

// Return a random element of Z/pZ.
int zp_rand();

// Compute base^exp mod p using modular exponentiation.
int zp_exp(int base, int exp);

// Return a n^th primitive root of unity in Z/pZ.
// We should have p = nq + 1 where n is a power of 2,
// and q is odd.
int zp_prim_root();

// Return the smallest n^th primitive root of unity in Z/pZ.
// We should have p = nq + 1 where n is a power of 2,
// and q is odd.
int zp_prim_root_min();

#endif