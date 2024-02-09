#ifndef zp_array_H
#define zp_array_H

#include <stdbool.h>

#include "finite_field.h"

typedef zp_t *zp_array;

// Some operations on zp_arrays of int.

/******************************************************/

/* Printing functions */

// Print the int of tab in stdout.
void zp_array_print(zp_array tab, int tab_size);

/******************************************************/

/* Initialization functions */

// Allocate memory of an zp_array of tab_size int.
zp_array zp_array_new(int tab_size);

// Same as above but set the zp_array elements to 0.
zp_array zp_array_new_zeros(int tab_size);

// Combine zp_array_new and zp_array_set.
zp_array zp_array_new_set(int tab_size, ...);

// Same as zp_array_new but fill it with random zp_t.
zp_array zp_array_new_rand(int tab_size);

/******************************************************/

/* Cleaning functions */

// Liberate memory occupied by tab.
void zp_array_free(zp_array tab);

/******************************************************/

/* Utility fonctions */

// Set an already allocate zp_array with the given values.
void zp_array_set(zp_array tab, int tab_size, ...);

// Compare two zp_array with the same lengh. Return true if
// they have the same elements, false else.
bool zp_array_equal(zp_array tab1, zp_array tab2, int tab_size);

// Add random errors on a given zp_array. Add at most
// number_errors errors
void zp_array_add_errors(zp_array tab, int tab_size, int number_errors);

/******************************************************/

/* FFT functions */

// Compute the fast fourier transform
// d is a power of 2
// coeffs should have < d elements
zp_array zp_array_fft(zp_array coeffs, int d);

// Compute the inverse of the fast fourrier transform
// d is a power of 2
// eval should have < d elements
zp_array zp_array_inv_fft(zp_array eval, int d);

/******************************************************/

/* Formal power series functions */

// Multiplication of formal series P and Q.
zp_array formal_serie_mul(zp_array P, zp_array Q, int d);

// Compute the d first coefficients of the inverse serie of P
// P should be at least of size d.
// d should be a power of 2 and 2*d should divide p-1.
zp_array formal_serie_inv(zp_array P, int d);

// deg_p is the degree of op_p, op_p should be at least of size deg_p + 1.
// deg_d is the degree of op_d, op_d should be at least of size deg_d + 1.
// It compute the quotient and rest in the euclidian division of op_p by op_d.
// It store them in quo and rem
// quo will be of degree deg_p - deg_d and of size deg_p - deg_d + 1.
// rem will be of degree at most deg_d - 1 and is always of size deg_d, padding with 0.
void formal_serie_fast_euc_div(zp_array *quo, zp_array *rem, zp_array op_p, int deg_p, zp_array op_d, int deg_d);

#endif