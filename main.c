#include <stdlib.h>
#include <stdio.h>

#include "rs_code.h"
#include "polynomial.h"
#include "finite_field.h"

void print_tab(int *T, int n)
{
	for (int i = 0; i < n; i++)
		printf("%d ", T[i]);
	printf("\n");
}

int main(int argc, char **argv)
{
	printf("Test des opération de base:\n");

	printf("\nLes polynômes:\n");
	poly *f = new_poly_from_coeffs(2, 5, 4, 3);
	print_poly(f);
	poly *g = new_poly_from_coeffs(2, 1, 3, 2);
	print_poly(g);

	printf("La somme:\n");
	poly *sum = new_poly();
	add_poly(sum, f, g);
	print_poly(sum);
	free_full_poly(sum);

	printf("La différence:\n");
	poly *dif = new_poly();
	sub_poly(dif, f, g);
	print_poly(dif);
	free_full_poly(dif);

	printf("Le produit:\n");
	poly *prod = new_poly();
	mul_poly(prod, f, g);
	print_poly(prod);
	free_full_poly(prod);

	printf("La multiplication par un scalaire:\n");
	int x = 10;
	poly *xf = new_poly();
	mul_scalar_poly(xf, x, f);
	print_poly(xf);
	free_full_poly(xf);

	printf("Le quotient et le reste:\n");
	poly *quo = new_poly();
	poly *rem = new_poly();
	euclid_div_poly(quo, rem, f, g);
	print_poly(quo);
	print_poly(rem);
	free_full_poly(quo);
	free_full_poly(rem);

	printf("Euclide étendu:\n");
	poly *d, *u, *v;
	xgcd_poly(&d, &u, &v, f, g);
	print_poly(d);
	print_poly(u);
	print_poly(v);
	free_full_poly(d);
	free_full_poly(u);
	free_full_poly(v);

	free_full_poly(f);
	free_full_poly(g);

	printf("\nTest de l'interpolation:\n");

	int n = 11;
	int a[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	f = new_poly_from_coeffs(n - 1, 231, 142, 229, 200, 0, 12, 2, 3, 123, 8, 28);
	print_poly(f);
	int b[n + 1];
	multi_eval_poly(f, n, a, b);
	free_full_poly(f);
	g = interpolation(a, b, n);
	print_poly(g);
	free_full_poly(g);

	printf("\nTest des codes de Reed-Solomon\n");
	n = 7;
	int A[] = {0, 1, 2, 3, 4, 5, 6};
	printf("Pour le message:\n");
	int k = 3;
	int M[] = {1, 2, 3};
	print_tab(M, k);
	printf("Le mot de code est:\n");
	int C[n];
	encoding(C, A, n, M, k);
	print_tab(C, n);
	printf("Le message reçu avec erreurs est:\n");
	int B[] = {1, 6, 123, 456, 57, 86, 121};
	print_tab(B, n);
	printf("Avec la méthode de Gao, on trouve:\n");
	int m[k];
	decoding(m, A, B, n, k);
	print_tab(m, k);
	return EXIT_SUCCESS;
}