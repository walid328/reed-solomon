#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

// Basic operations in Z/pZ.
// Elements of Z/pZ are represented by 0, 1, ..., p-1.

typedef unsigned int zp_t;

// Compute a mod p.
zp_t zp_mod(int a);

// Compute a + b mod p.
zp_t zp_add(zp_t a, zp_t b);

// Compute a - b mod p.
zp_t zp_sub(zp_t a, zp_t b);

// Compute -a mod p.
zp_t zp_opp(zp_t a);

// Compute a*b mod p.
zp_t zp_mul(zp_t a, zp_t b);

// Compute a^(-1) mod p via xgcd.
zp_t zp_inv_xgcd(zp_t a);

// Return a^(-1) mod p via precomputed table.
zp_t zp_inv(zp_t a);

// Return a random element of Z/pZ.
zp_t zp_rand();

// Compute base^exp mod p using modular exponentiation.
zp_t zp_exp(zp_t base, zp_t exp);

// Return a n^th primitive root of unity in Z/pZ.
// We should have p = nq + 1 where n is a power of 2,
// and q is odd.
zp_t zp_prim_root();

// Return the smallest n^th primitive root of unity in Z/pZ.
// We should have p = nq + 1 where n is a power of 2,
// and q is odd.
zp_t zp_prim_root_min();

#endif