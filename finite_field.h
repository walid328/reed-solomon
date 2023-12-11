#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

struct
{
    int degree;
    int *coefficients;
} typedef poly;

int size_int_str(int n);

void str_add_int(char* string, int *index, int n);

char *string_polynomial_iso_length(poly *q);

char *string_polynomial_minimal(poly *q);

void print_polynomial_iso_length(poly *q);

void print_polynomial_minimal(poly *q);

void print_polynomial(poly *q);

poly *add_polynomials(poly *p1, poly *p2);

poly *mul_polynomials(poly *p1, poly *p2);

void euclid_division(poly *p1, poly *p2, poly *q, poly *r);

poly *xgcd(poly *p1, poly *p2, poly *u, poly *v);

poly *interpolation(int *a, int *b);

#endif