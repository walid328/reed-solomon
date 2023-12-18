#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

typedef struct polynomial poly;

int size_int_str(int n);

void str_add_int(char *string, int *index, int n);

char *str_poly_iso_length(poly *q);

char *str_poly_min(poly *q);

void print_poly_iso_length(poly *q);

void print_poly_min(poly *q);

void print_poly(poly *q);

poly *new_poly(void);

poly *new_poly_from_coeffs(int deg, ...);

poly *new_poly_0(void);

poly *new_poly_1(void);

poly *new_poly_from_copy(poly *source);

void set_poly(poly *q, int degree, int *coefficients);

void free_poly(poly *q);

void free_poly_full(poly *q);

poly *add_poly(poly *p1, poly *p2);

poly *substract_poly(poly *p1, poly *p2);

poly *mul_poly(poly *p1, poly *p2);

poly *mul_poly_scalar(poly *q, int n);

poly *derivate(poly *q);

int evaluate(poly *q, int x);

int *multi_evaluate(poly *q, int *a, int n);

// return the inverse of n in Z/pZ
int inverse_zp(int n);

void euclid_division(poly *p1, poly *p2, poly *q, poly *r);

poly *xgcd(poly *p1, poly *p2, poly *u, poly *v);

poly *interpolation(int *a, int *b, int n);

#endif