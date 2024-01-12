#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

// Basic implementation of Fp[x].

typedef struct polynomial poly;

/******************************************************/

void poly_set(poly *f, int degree, int *coefficients);

int poly_degree(poly *f);

int *poly_coeffs(poly *f);

/******************************************************/

/* Printing functions */

int size_int_str(int n);

void str_add_int(char *string, int *index, int n);

char *str_poly_iso_length(poly *q);

void poly_print(poly *q);

/******************************************************/

/* Initialization functions */

// Returns a zero polynomial.
// The zero polynomial is represented with degree = -1
// and a coeffs = NULL.
poly *poly_new(void);

// Returns a polynomial of degree deg,
// the coefficients are given in arguments.
// For instance for "x^2 + x + 1" the call would
// be poly_new_from_poly_coeffs(2, 1, 1, 1).
poly *poly_new_from_poly_coeffs(int deg, ...);

// Returns the polynomial equal to 1.
poly *poly_new_1(void);

int next_coeff(char *str, int *i);

int *coeffs_from_str(char *str, int *deg);

// Returns a polynomial from a string given
// in argument. For instance: "1 + x + x^2".
poly *poly_new_from_str(char *str);

// Returns a polynomial equal to a polynomial
// given in argument.
poly *poly_new_from_copy(poly *source);

// Returns a random polynomial.
poly *poly_new_rand(int deg);

/******************************************************/

/* Cleaning functions */

// Resets the polynomial f.
void poly_clear(poly *f);

// Conveniance function to resets the coeffs too.
void poly_clear_full(poly *f);

// Free a polynomial.
void poly_free(poly *f);

// Conveniance function to free the coeffs too.
void poly_free_full(poly *f);

/******************************************************/

/* Basic operations */

// Adds op1 and op2 and stores the sum in rop.
void poly_add(poly *rop, poly *op1, poly *op2);

// Substracts op1 and op2 and stores the difference in rop.
void poly_sub(poly *rop, poly *op1, poly *op2);

// Multiplies op1 and op2 and stores the product in rop.
void poly_mul(poly *rop, poly *op1, poly *op2);

// Multiplies op1 and op2 and stores the product in rop.
void poly_mul_scalar(poly *rop, int op1, poly *op2);

// Computes the euclidian division of op1 by op2, i.e,
// op1 = q*op2 + r where deg(r) < deg(op2).
void poly_euc_div(poly *q, poly *r, poly *op1, poly *op2);

// Extended Euclidean algorithm.
// d = u*op1 + v*op2.
// d, u, v shouldn't be initialized.
// d is a monic polynomial.
void poly_xgcd(poly **d, poly **u, poly **v, poly *op1, poly *op2);

// Extended Euclidean algorithm but stops when remainder has a degree smaller
// than n.
// d = u*op1 + v*op2.
// d, u, v shouldn't be initialized.
void poly_xgcd_partial(poly **d, poly **u, poly **v, poly *op1, poly *op2, int n);

/******************************************************/

/* DFT functions */

// Derives op and stores the derivative in rop.
void poly_deriv(poly *rop, poly *op);

// Returns y = f(x).
int poly_eval(poly *f, int x);

// Evaluates the polynomial f for n value in the
// tab a and stores the result in tab b.
void poly_eval_multi(poly *f, int n, int *a, int *b);

// dft for polynomials
void poly_dft(poly *f, int **eval);

// Computes a polynomial using n points (a, b) such
// that f(a) = b. We use Lagrange polynomials to
// compute the inverse of the Vandermonde matrix.
poly *interpolation(int *a, int *b, int n);

/******************************************************/

/* FFT functions */

// Works if p = q*d + 1 where d is a power of 2

// "Splits" the array tab of size tab_size into two arrays
// containing elements of even indexes and even indexes.
void array_split_all(int **even, int *even_size, int **odd, int *odd_size, int *tab, int tab_size);

// splits array with even length
void array_split(int **even, int **odd, int *tab, int tab_size);

// inverse of split array
void array_merge(int **tab, int *even, int *odd, int subtab_size);

// "Splits" the polynomial f in two polynomial even and odd
// such that f(x) = even(x^2) + x*odd(x^2).
// f should have a degree of d-1 at most.
void poly_split(poly *even, poly *odd, poly *f);

// fast fourier transform for d a power of twoo
// omega should be a d^th primitive root of unity
int *fft(int *f, int d, int omega);

// fft for polynomials
void poly_fft(poly *f, int **eval);

// inverse fft in int* format
int *inv_fft(int *eval, int d, int omega);

// inverse fft in poly format
void poly_inv_fft(poly *f, int *eval);

// Polynomial multiplication using fft.
void poly_mul_fft(poly *rop, poly *op1, poly *op2);

#endif