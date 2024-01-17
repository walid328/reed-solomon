#include <stdlib.h>
#include <stdio.h>

#include "rs_code.h"
#include "polynomial.h"
#include "finite_field.h"
#include "field_spec.h"

int p = 193;
int q = 0;
int n = 0;
int omega = 0;
int *omegas = NULL;

void array_print(int *tab, int tab_size)
{
	for (int i = 0; i < tab_size; i++)
		printf("%d ", tab[i]);
	printf("\n");
}

int main(int argc, char **argv)
{
	update_field_spec();
	printf("\nTest des codes de Reed-Solomon\n");
	int eval_size = 7;
	int A[] = {0, 1, 2, 3, 4, 5, 6};
	printf("Pour le message:\n");
	int message_size = 3;
	int M[] = {1, 2, 3};
	array_print(M, message_size);
	printf("Le mot de code est:\n");
	int C[eval_size];
	rs_encode(C, A, eval_size, M, message_size);
	array_print(C, eval_size);
	printf("Le message reçu avec erreurs est:\n");
	int B[] = {1, 6, 123, 456, 57, 86, 121};
	array_print(B, eval_size);
	printf("Avec la méthode de Gao, on trouve:\n");
	int m[message_size];
	rs_decode(m, A, B, eval_size, message_size);
	array_print(m, message_size);

	printf("\nTest de l'exponentiation modulaire:\n");
	int base = 146;
	int exp = 151;
	int res = zp_exp(base, exp);
	printf("%d^%d mod %d = %d\n", base, exp, p, res);

	printf("\nTest de racine primitive de l'unité\n");
	int G = zp_prim_root();
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
	poly *t = poly_new_from_coeffs(10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
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
	printf("q : %d , d : %d\n", q, n);
	printf("omega : %d\n", omega);
	poly *fftest = poly_new_rand(10);
	int *eval_fft;
	int *eval_dft;
	poly_fft(fftest, &eval_fft);
	poly_dft(fftest, &eval_dft);
	array_print(eval_fft, n);
	array_print(eval_dft, n);

	poly *test_fftest = poly_new();
	poly_inv_fft(test_fftest, eval_fft);
	poly_print(fftest);
	poly_print(test_fftest);

	poly_free_full(fftest);
	poly_free_full(test_fftest);
	free(eval_fft);
	free(eval_dft);

	printf("\nTest serie formelle:\n");
	int d = 64;
	poly *serie_A = poly_new_rand(n - 1);
	poly *serie_B = poly_new_rand(n - 1);

	int *inv_serie_A = formal_serie_inv(poly_coeffs(serie_A), d);
	array_print(poly_coeffs(serie_A), d / 2);
	int *serie = formal_serie_inv(inv_serie_A, d);
	array_print(serie, d / 2);
	free(inv_serie_A);
	free(serie);

	poly_free_full(serie_A);
	poly_free_full(serie_B);

	free(omegas);
	return EXIT_SUCCESS;
}