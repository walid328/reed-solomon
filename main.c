#include <stdlib.h>
#include <stdio.h>

#include "finite_field.h"

int main(int argc, char **argv)
{
	int coeffs_p1[] = {5, 4, 3};
	int coeffs_p2[] = {1, 3, 2};

	poly *p1 = new_poly();
	set_poly(p1, 2, coeffs_p1);
	poly *p2 = new_poly();
	set_poly(p2, 2, coeffs_p2);

	poly *sum = add_poly(p1, p2);
	print_poly(sum);

	poly *prod = mul_poly(p1, p2);
	print_poly(prod);

	free_poly(p1);
	free_poly(p2);
	free_poly_full(sum);
	free_poly_full(prod);

	int coeffs_p3[] = {1, 3, 2, 19};
	poly *p3 = new_poly();
	set_poly(p3, 3, coeffs_p3);

	int coeffs_p4[] = {5, 2};
	poly *p4 = new_poly();
	set_poly(p4, 1, coeffs_p4);

	poly *quo = new_poly();
	poly *rest = new_poly();
	euclid_division(p3, p4, quo, rest);
	print_poly_iso_length(p3);
	print_poly_iso_length(p4);
	print_poly_iso_length(quo);
	print_poly_iso_length(rest);
	free_poly(p3);
	free_poly(p4);
	free_poly_full(quo);
	free_poly_full(rest);

	int n = 10;
	int a[] = {200, 100, 49, 1, 256, 0, 67, 98, 212, 300, 27};
	int coeffs_q[] = {52, 298, 189, 100, 167, 234, 232, 10, 0, 1, 56};
	poly *q = new_poly();
	set_poly(q, 10, coeffs_q);
	int *b = multi_evaluate(q, a, n);
	print_poly_iso_length(q);
	poly *qq = interpolation(a, b, n);
	print_poly_iso_length(qq);
	free_poly(q);
	free_poly_full(qq);
	free(b);

	poly *p5 = new_poly_from_str(" 57 + 193*x +  83*x^2 + 197*x^3 + 200*x^4 + 124*x^5 +  16*x^6 + 246*x^7 + 137*x^8 + 300*x^9 + 167*x^10");
	print_poly_iso_length(p5);
	free_poly_full(p5);

	poly *p6 = new_rand_poly(5);
	print_poly_iso_length(p6);
	free_poly_full(p6);

	return EXIT_SUCCESS;
}