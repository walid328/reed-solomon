#include <stdlib.h>
#include <stdio.h>

#include "rs_code.h"
#include "polynomial.h"
#include "finite_field.h"

void array_print(int *tab, int tab_size)
{
	for (int i = 0; i < tab_size; i++)
		printf("%d ", tab[i]);
	printf("\n");
}

int main(int argc, char **argv)
{
	int Q, D;
	for (Q = p - 1, D = 1; (Q & 1) == 0; Q >>= 1, D <<= 1)
		;
	printf("Test des opération de base:\n");

	printf("\nLes polynômes:\n");
	poly *f = poly_new_from_poly_coeffs(2, 5, 4, 3);
	poly_print(f);
	poly *g = poly_new_from_poly_coeffs(2, 1, 3, 2);
	poly_print(g);

	printf("La somme:\n");
	poly *sum = poly_new();
	poly_add(sum, f, g);
	poly_print(sum);
	poly_free_full(sum);

	printf("La différence:\n");
	poly *dif = poly_new();
	poly_sub(dif, f, g);
	poly_print(dif);
	poly_free_full(dif);

	printf("Le produit:\n");
	poly *prod = poly_new();
	poly *prod_fft = poly_new();
	poly_mul(prod, f, g);
	poly_mul_fft(prod_fft, f, g);
	poly_print(prod);
	poly_print(prod_fft);
	poly_free_full(prod);
	poly_free_full(prod_fft);

	printf("La multiplication par un scalaire:\n");
	int x = 10;
	poly *xf = poly_new();
	poly_mul_scalar(xf, x, f);
	poly_print(xf);
	poly_free_full(xf);

	printf("Le quotient et le reste:\n");
	poly *quo = poly_new();
	poly *rem = poly_new();
	poly_euc_div(quo, rem, f, g);
	poly_print(quo);
	poly_print(rem);
	poly_free_full(quo);
	poly_free_full(rem);

	printf("Euclide étendu:\n");
	poly *d, *u, *v;
	poly_xgcd(&d, &u, &v, f, g);
	poly_print(d);
	poly_print(u);
	poly_print(v);
	poly_free_full(d);
	poly_free_full(u);
	poly_free_full(v);

	poly_free_full(f);
	poly_free_full(g);

	printf("\nTest de l'interpolation:\n");

	int n = 11;
	int a[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	f = poly_new_from_poly_coeffs(n - 1, 231, 142, 229, 200, 0, 12, 2, 3, 123, 8, 28);
	poly_print(f);
	int b[n + 1];
	poly_eval_multi(f, n, a, b);
	poly_free_full(f);
	g = interpolation(a, b, n);
	poly_print(g);
	poly_free_full(g);

	printf("\nTest des codes de Reed-Solomon\n");
	n = 7;
	int A[] = {0, 1, 2, 3, 4, 5, 6};
	printf("Pour le message:\n");
	int k = 3;
	int M[] = {1, 2, 3};
	array_print(M, k);
	printf("Le mot de code est:\n");
	int C[n];
	rs_encode(C, A, n, M, k);
	array_print(C, n);
	printf("Le message reçu avec erreurs est:\n");
	int B[] = {1, 6, 123, 456, 57, 86, 121};
	array_print(B, n);
	printf("Avec la méthode de Gao, on trouve:\n");
	int m[k];
	rs_decode(m, A, B, n, k);
	array_print(m, k);

	printf("\nTest de l'exponentiation modulaire:\n");
	int base = 146;
	int exp = 151;
	int res = zp_exp(base, exp);
	printf("%d^%d mod %d = %d\n", base, exp, p, res);

	printf("\nTest de racine primitive de l'unité\n");
	int G = zp_prim_root(Q, D);
	printf("%d\n", G);

	printf("\nTest split array:\n");
	int tab_size = 19;
	int *tab = malloc(tab_size * sizeof(int));
	for (int i = 0; i < tab_size; i++)
		tab[i] = i;
	array_print(tab, tab_size);
	int *even_tab = NULL;
	int even_size = 0;
	int *odd_tab = NULL;
	int odd_size = 0;
	array_split_all(&even_tab, &even_size, &odd_tab, &odd_size, tab, tab_size);
	array_print(even_tab, even_size);
	array_print(odd_tab, odd_size);
	free(even_tab);
	free(odd_tab);
	free(tab);

	printf("\nTest split poly:\n");
	poly *t = poly_new_from_poly_coeffs(10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
	poly_print(t);
	poly *even = poly_new();
	poly *odd = poly_new();
	poly_split(even, odd, t);
	poly_print(even);
	poly_print(odd);
	poly_free_full(t);
	poly_free_full(even);
	poly_free_full(odd);

	printf("\nTest fft:\n");
	int omega = zp_prim_root_min(Q, D);
	printf("q : %d , d : %d\n", Q, D);
	printf("omega : %d\n", omega);
	poly *fftest = poly_new_rand(10);
	int *eval_fft;
	int *eval_dft;
	poly_fft(fftest, &eval_fft);
	poly_dft(fftest, &eval_dft);
	array_print(eval_fft, D);
	array_print(eval_dft, D);

	poly *test_fftest = poly_new();
	poly_inv_fft(test_fftest, eval_fft);
	poly_print(fftest);
	poly_print(test_fftest);

	poly_free_full(fftest);
	poly_free_full(test_fftest);
	free(eval_fft);
	free(eval_dft);

	return EXIT_SUCCESS;
}