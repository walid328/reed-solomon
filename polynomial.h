#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

// Basic implementation of Fp[x].

typedef struct polynomial poly;

/******************************************************/

void set_poly(poly *q, int degree, int *coefficients);

/******************************************************/

/* Printing functions */

int size_int_str(int n);

void str_add_int(char *string, int *index, int n);

char *str_poly_iso_length(poly *q);

char *str_poly_min(poly *q);

void print_poly_iso_length(poly *q);

void print_poly_min(poly *q);

void print_poly(poly *q);

/******************************************************/

/* Initialization functions */

// Returns an empty polynomial.
poly *new_poly(void);

// Returns a polynomial of degree deg,
// the coefficients are given in arguments.
// For instance for "x^2 + x + 1" the call would
// be new_poly_from_coeffs(2, 1, 1, 1).
poly *new_poly_from_coeffs(int deg, ...);

// Return the zero polynomial.
poly *new_poly_0(void);

// Returns the polynomial equal to 1.
poly *new_poly_1(void);

int next_coeff(char *str, int *i);

int *coeffs_from_str(char *str, int *deg);

// Returns a polynomial from a string given
// in argument. For instance: "1 + x + x^2".
poly *new_poly_from_str(char *str);

// Returns a polynomial equal to a polynomial
// given in argument.
poly *new_poly_from_copy(poly *source);

// Returns a random polynomial.
poly *new_rand_poly(int deg);

/******************************************************/

/* Cleaning functions */

// Resets the polynomial f.
void clear_poly(poly *f);

// Conveniance function to resets the coeffs too.
void clear_full_poly(poly *f);

// Free a polynomial.
void free_poly(poly *f);

// Conveniance function to free the coeffs too.
void free_full_poly(poly *f);

/******************************************************/

/* Basic operations */

// Adds op1 and op2 and stores the sum in rop.
void add_poly(poly *rop, poly *op1, poly *op2);

// Substracts op1 and op2 and stores the difference in rop.
void sub_poly(poly *rop, poly *op1, poly *op2);

// Multiplies op1 and op2 and stores the product in rop.
void mul_poly(poly *rop, poly *op1, poly *op2);

// Multiplies op1 and op2 and stores the product in rop.
void mul_scalar_poly(poly *rop, int op1, poly *op2);

// Computes the euclidian division of op1 by op2, i.e,
// op1 = q*op2 + r where deg(r) < deg(op2).
void euclid_div_poly(poly *q, poly *r, poly *op1, poly *op2);

// Extended Euclidean algorithm.
// d = u*op1 + v*op2.
// d, u, v shouldn't be initialized.
// d is a monic polynomial.
void xgcd_poly(poly **d, poly **u, poly **v, poly *op1, poly *op2);

/******************************************************/

/* DFT functions */

// Derives op and stores the derivative in rop.
void deriv_poly(poly *rop, poly *op);

// Returns y = f(x).
int eval_poly(poly *f, int x);

// Evaluates the polynomial f for n value in the
// tab a and stores the result in tab b.
void multi_eval_poly(poly *f, int n, int *a, int *b);

// Computes a polynomial using n points (a, b) such
// that f(a) = b. We use Lagrange polynomials to
// compute the inverse of the Vandermonde matrix.
poly *interpolation(int *a, int *b, int n);

#endif