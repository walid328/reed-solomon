#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <stdbool.h>

#include "array.h"

// Basic implementation of Fp[x].

typedef struct polynomial *poly;

/******************************************************/

/* Access functions */

// Return deg(f).
int poly_deg(const poly f);

// Return the array of coefficients of f.
array poly_coeffs(const poly f);

// Set the polynomial f with the given values.
// coeffs should be manually allocated (with
// array_new for instance).
void poly_set(poly f, int deg, array coeffs);

/******************************************************/

/* Printing functions */

// Print f in stdout.
void poly_print(const poly f);

/******************************************************/

/* Initialization functions */

// Return a zero polynomial, i.e with degree -1 and
// NULL coefficient.
poly poly_new(void);

// Return a polynomial with degree deg, the coefficients
// are given in arg. For instance poly_new_set(2, 1, 1, 1)
// returns "x^2 + x + 1"
poly poly_new_set(int deg, ...);

// Return a polynomial of degree deg.
// Allocate the memory for the polynomial and for his coefficients
// Doesn't set the coefficients.
poly poly_new_deg(int deg);

// Return a polynomial from a string in arg.
poly poly_new_str(char *str);

// Return a random polynomial with degree deg.
poly poly_new_rand(int deg);

/******************************************************/

/* Cleaning functions */

// Reset the polynomial f to zero.
void poly_clear(poly f);

// Reset the qty polynomials passed as arguments.
void poly_clear_multi(int qty, ...);

// Liberate the memory occupied by f.
void poly_free(poly f);

// Liberate the memory occupied by the qty polynomials
// passed as arguments.
void poly_free_multi(int qty, ...);

/******************************************************/

/* Utility functions */

// Return true if polynomials are the same, else false.
bool poly_equal(const poly f, const poly g);

// Put a copy of the source polynomial in destination.
void poly_copy(poly dst, const poly src);

// Put the reverse of a given polynomial in rop.
// For instance with "x^3 + 4x^2 + 7" will
// return "7x^3 + 4x + 1".
void poly_rev(poly rop, const poly f);

// Compute the derivate polynomial of op and store it in rop.
void poly_deriv(poly rop, const poly op);

/******************************************************/

/* Basic operations */

// Compute op1 + op2 and store it in rop.
void poly_add(poly rop, const poly op1, const poly op2);

// Compute op1 - op2 and store it in rop.
void poly_sub(poly rop, const poly op1, const poly op2);

// Compute op1 * op2 and store it in rop.
void poly_mul(poly rop, const poly op1, const poly op2);

// Compute op1 * op2 and store it in rop.
void poly_mul_scalar(poly rop, int op1, const poly op2);

// Compute the quotient and remainder of the euclidian
// division of op1 and op2 and store them in q and r.
// op1 = q * op2 + r, deg(r) < deg(op1)
void poly_euc_div(poly q, poly r, const poly op1, const poly op2);

// Compute the extended euclidian algorithm with op1 and op2.
// u*op1 + v*op2 = d
void poly_xgcd(poly *d, poly *u, poly *v, const poly op1, const poly op2);

// Same as above but stop when deg(d) < limit.
void poly_xgcd_partial(poly *d, poly *u, poly *v, const poly op1, const poly op2, int limit);

/******************************************************/

/*DFT functions */

// Compute f(a).
int poly_eval(const poly f, int a);

// Compute f(a) for each a in tab.
array poly_eval_array(const poly f, const array tab, int tab_size);

// Same as above but with powers of a privitime root of unity.
// f should have a degree < d.
array poly_dft(const poly f);

// Find the polynomial f such that f(a_i) = b_i for a in points
// and b in eval using lagrange polynomials.
void interpolation(poly rop, const array points, const array eval, int size);

// Same as above but the evalutation points are the powers of a
// primitive root of unity.
void poly_inv_dft(poly rop, array eval);

/******************************************************/

/* FFT functions */

// Works if p = q*d + 1 where d is a power of 2

// Fast fourrier transform for polynomials. Same as poly_dft but
// faster.
array poly_fft(const poly f);

// Inverse of fft. Same as interpolation but faster.
void poly_inv_fft(poly rop, array eval);

/******************************************************/

/* Fast operations */

// Works if p = q*d + 1 where d is a power of 2

// Same as poly_mul but faster.
void poly_fast_mul(poly rop, const poly op1, const poly op2);

// Same as poly_euc_div but faster.
void poly_fast_euc_div(poly q, poly r, const poly op1, const poly op2);

poly *poly_half_gcd(const poly r0, const poly r1);

poly *poly_fast_gcd_matrix(const poly r0, const poly r1);

void poly_fast_xgcd(poly d, poly u, poly v, const poly op_a, const poly op_b);

#endif