#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <stdbool.h>

#include "zp_array.h"
#include "finite_field.h"

// Basic implementation of Fp[x].

typedef struct polynomial *poly;

/******************************************************/

/* Access functions */

// Return deg(f).
int poly_deg(const poly f);

// Return the zp_array of coefficients of f.
zp_array poly_coeffs(const poly f);

// Return the leading coefficient of f.
zp_t poly_leading_coeff(const poly f);

// Set the polynomial f with the given values.
// coeffs should be manually allocated (with
// zp_array_new for instance).
void poly_set(poly f, int deg, zp_array coeffs);

/******************************************************/

/* Printing functions */

// Print f in stdout.
void poly_print(const poly f);

/******************************************************/

/* Initialization functions */

// Return a zero polynomial, i.e with degree -1 and
// NULL coefficient.
poly poly_new(void);

// Return a polynomial of degree deg.
// Allocate the memory for the polynomial and for his coefficients
// Doesn't set the coefficients.
poly poly_new_deg(int deg);

// Return a polynomial with degree deg, the coefficients
// are given in arg. For instance poly_new_set(2, 3, 2, 1)
// returns "x^2 + 2x + 3"
poly poly_new_set(int deg, ...);

// Return a polynomial from a string in arg.
poly poly_new_str(char *str);

// Return a random polynomial with degree deg.
poly poly_new_rand(int deg);

// Return a copy of the source polynomial.
poly poly_new_copy(const poly src);

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

// Return true if f = 0, else false.
bool poly_is_zero(const poly f);

// Return true if polynomials are the same, else false.
bool poly_equal(const poly f, const poly g);

// Put a copy of the source polynomial in destination.
void poly_copy(poly dst, const poly src);

// Clear f, set his degree and allocate his zp_array coeffs.
void poly_reset_deg(poly f, int deg);

// Clear f and set it like in poly_new_set.
void poly_reset_coeffs(poly f, int deg, ...);

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
void poly_mul_scalar(poly rop, zp_t op1, const poly op2);

// Compute the quotient and remainder of the euclidian
// division of op1 and op2 and store them in q and r.
// op1 = q * op2 + r, deg(r) < deg(op1)
void poly_euc_div(poly q, poly r, const poly op1, const poly op2);

// Compute the extended euclidian algorithm with op1 and op2.
// u*op1 + v*op2 = d
void poly_xgcd(poly d, poly u, poly v, const poly op1, const poly op2);

// Same as above but stop when deg(d) < limit.
void poly_xgcd_partial(poly d, poly u, poly v, const poly op1, const poly op2, int limit);

/******************************************************/

/*DFT functions */

// Compute f(a).
int poly_eval(const poly f, int a);

// Compute f(a) for each a in tab.
zp_array poly_eval_zp_array(const poly f, const zp_array tab, int tab_size);

// Find the polynomial f such that f(a_i) = b_i for a in points
// and b in eval using lagrange polynomials.
void interpolation(poly rop, const zp_array points, const zp_array eval, int size);

// Discrete Fourrier transform for polynomials.
// The evaluation points are powers of a primitive d-th root of unity.
// d is a power of 2 dividing p-1.
zp_array poly_dft(const poly f, int d);

// Inverse of fft.
void poly_inv_dft(poly rop, zp_array eval, int d);

/******************************************************/

/* FFT functions */

// Works if p = q*d + 1 where d is a power of 2

// Fast Fourrier transform for polynomials. Same as poly_dft but in O(n log(n)).
// The evaluation points are powers of a primitive d-th root of unity.
// d is a power of 2 dividing p-1.
// d is at least deg(f) + 1.
zp_array poly_fft(const poly f, int d);

// Inverse of fft. Same as inv_dft with a better time complexity.
void poly_inv_fft(poly rop, zp_array eval, int d);

/******************************************************/

/* Fast operations */

// Works if p = q*d + 1 where d is a power of 2

// Same as poly_mul with a better time complexity.
void poly_fast_mul(poly rop, const poly op1, const poly op2);

// Same as poly_euc_div with a better time complexity.
void poly_fast_euc_div(poly q, poly r, const poly op1, const poly op2);

// Same as poly_xgcd with a better time complexity.
void poly_fast_xgcd(poly d, poly u, poly v, const poly op1, const poly op2);

// Same as poly_xgcd_partial with a better time complexity.
void poly_fast_xgcd_partial(poly d, poly u, poly v, const poly op1, const poly op2, int limit);

#endif