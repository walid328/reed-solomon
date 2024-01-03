#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H
#define p 193

// Basic operations in Z/pZ.

// Returns the representant of n in Z/pZ.
int mod_zp(int n);

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

// Returns x^n mod p using modular exponentiation.
int exp_zp(int base, int exp);

// Returns a d^th primitive root of unity in Z/pZ.
// We should have p = dq + 1 where d is a power of 2,
// and q is odd.
int primitive_root_zp(int q, int d);

#endif