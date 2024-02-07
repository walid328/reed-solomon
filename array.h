#ifndef ARRAY_H
#define ARRAY_H

#include <stdbool.h>

#include "finite_field.h"

typedef zp_t *array;

// Some operations on arrays of int.

/******************************************************/

/* Printing functions */

// Print the int of tab in stdout.
void array_print(array tab, int tab_size);

/******************************************************/

/* Initialization functions */

// Allocate memory of an array of tab_size int.
array array_new(int tab_size);

// Same as above but set the array elements to 0.
array array_new_zeros(int tab_size);

// Combine array_new and array_set.
array array_new_set(int tab_size, ...);

// Same as array_new but fill it with random zp_t.
array array_new_rand(int tab_size);

/******************************************************/

/* Cleaning functions */

// Liberate memory occupied by tab.
void array_free(array tab);

/******************************************************/

/* Utility fonctions */

// Set an already allocate array with the given values.
void array_set(array tab, int tab_size, ...);

// Compare two array with the same lengh. Return true if
// they have the same elements, false else.
bool array_equal(array tab1, array tab2, int tab_size);

// Add random errors on a given array.
void array_add_errors(array tab, int tab_size, int number_errors);

/******************************************************/

/* FFT functions */

// Compute the fast fourier transform
// d is a power of 2
// coeffs should have < d elements
array array_fft(array coeffs, int d);

// Compute the inverse of the fast fourrier transform
// d is a power of 2
// eval should have < d elements
array array_inv_fft(array eval, int d);

/******************************************************/

/* Formal power series functions */

// Multiplication of formal series P and Q.
array formal_serie_mul(array P, array Q, int d);

// Compute the d first coefficients of the inverse serie of P
// P should be at least of size d.
// d should be a power of 2 and 2*d should divide p-1.
array formal_serie_inv(array P, int d);

// deg_p is the degree of op_p, op_p should be at least of size deg_p + 1.
// deg_d is the degree of op_d, op_d should be at least of size deg_d + 1.
// It compute the quotient and rest in the euclidian division of op_p by op_d.
// It store them in quo and rem
// quo will be of degree deg_p - deg_d and of size deg_p - deg_d + 1.
// rem will be of degree at most deg_d - 1 and is always of size deg_d, padding with 0.
void formal_serie_fast_euc_div(array *quo, array *rem, array op_p, int deg_p, array op_d, int deg_d);

#endif