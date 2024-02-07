#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "field_settings.h"
#include "finite_field.h"
#include "array.h"
#include "polynomial.h"
#include "rs_code.h"

int p = 0;
int q = 0;
int n = 0;
int omega = 0;
array omegas = NULL;
array inverses = NULL;

bool comp_dft(int qqt);

bool comp_inv_dft(int qqt);

bool comp_mul(int qqt);

bool comp_euc_div(int qqt);

bool comp_gcd(int qqt);

int main(void)
{
	field_settings_set(193);
	
	bool bon;
	bon = comp_dft(3);
	printf("%d comp_dft fait %d\n", bon, bon);
	bon = comp_inv_dft(3);
	printf("%d comp_inv_dft fait %d\n", bon, bon);
	bon = comp_mul(5);
	printf("%d comp_mul fait %d\n", bon, bon);
	bon = comp_euc_div(5);
	printf("%d comp_euc_div fait %d\n", bon, bon);
	bon = comp_gcd(5);
	printf("%d comp_gcd fait %d\n", bon, bon);
	
	field_settings_free();
	printf("fait\n");
	return EXIT_SUCCESS;
}

bool comp_dft(int qqt)
{
	bool bon = true;
	for (int dega = -1; dega < n; dega++)
	{
		for (int i = 0; i < qqt; i++)
		{
			poly pa = poly_new_rand(dega);
			int d = 1;
			while (d < dega + 1)
				d *= 2;
			while (d <= n)
			{
				array eval1 = poly_dft(pa, d);
				array eval2 = poly_fft(pa, d);
				if (!array_equal(eval1, eval2, d))
				{
					bon = false;
					printf("resultats distincts dft et fft\n");
					poly_print(pa);
					array_print(eval1, d);
					array_print(eval2, d);
				}
				array_free(eval1);
				array_free(eval2);
				d *= 2;
			}
			poly_free(pa);
		}
	}
	return bon;
}

bool comp_inv_dft(int qqt)
{
	bool bon = true;
	for (int d = 1; d <= n; d *= 2)
	{
		for (int i = 0; i < qqt; i++)
		{
			array eval = array_new_rand(d);
			poly pa = poly_new();
			poly pb = poly_new();
			poly_inv_dft(pa, eval, d);
			poly_inv_fft(pb, eval, d);
			if (!poly_equal(pa, pb))
			{
				bon = false;
				printf("resultats distincts inv_dft et inv_fft\n");
				array_print(eval, d);
				poly_print(pa);
				poly_print(pb);
			}
			array_free(eval);
			poly_free_multi(2, pa, pb);
		}
	}
	return bon;
}

bool comp_mul(int qqt)
{
	bool bon = true;
	for (int dega = -1; dega <= 2 * n; dega++)
	{
		for (int degb = -1; dega + degb + 1 <= 2 * n; degb++)
		{
			for (int i = 0; i < qqt; i++)
			{
				poly pa = poly_new_rand(dega);
				poly pb = poly_new_rand(degb);
				poly prod1 = poly_new();
				poly prod2 = poly_new();
				poly prod3 = poly_new();
				poly prod4 = poly_new();
				poly_mul(prod1, pa, pb);
				poly_fast_mul(prod2, pa, pb);
				poly_mul(prod3, pb, pa);
				poly_fast_mul(prod4, pb, pa);
				if (!poly_equal(prod1, prod2) || !poly_equal(prod1, prod3) || !poly_equal(prod1, prod4))
				{
					bon = false;
					printf("resultats distincts mul et fast_mul\n");
					printf("%d\n", poly_equal(prod1, prod2));
					printf("%d\n", poly_equal(prod1, prod3));
					printf("%d\n", poly_equal(prod1, prod4));
					poly_print(pa);
					poly_print(pb);
					poly_print(prod1);
					poly_print(prod2);
					poly_print(prod3);
					poly_print(prod4);
					fprintf(stderr, "1 is not inversible!\n");
					exit(EXIT_FAILURE);
				}
				poly_free_multi(6, pa, pb, prod1, prod2, prod3, prod4);
			}
		}
	}
	return bon;
}

bool comp_euc_div(int qqt)
{
	bool bon = true;
	for (int dega = -1; dega <= 2 * n - 1; dega++)
	{
		int min_degb = 0;
		if (dega + 1 - n > min_degb)
			min_degb = dega + 1 - n;
		for (int degb = min_degb; degb <= n; degb++)
		{
			for (int i = 0; i < qqt; i++)
			{
				poly pa = poly_new_rand(dega);
				poly pb = poly_new_rand(degb);
				poly quo1 = poly_new();
				poly rem1 = poly_new();
				poly quo2 = poly_new();
				poly rem2 = poly_new();
				poly_euc_div(quo1, rem1, pa, pb);
				poly_fast_euc_div(quo2, rem2, pa, pb);
				poly verif = poly_new();
				poly_mul(verif, quo1, pb);
				poly_add(verif, verif, rem1);
				if (!poly_equal(quo1, quo2) || !poly_equal(rem1, rem2) || !poly_equal(verif, pa) || poly_deg(rem1) >= poly_deg(pb))
				{
					bon = false;
					printf("resultats distincts euc_div et fast_euc_div\n");
					printf("pa :\n");
					poly_print(pa);
					printf("pb :\n");
					poly_print(pb);
					printf("quo1 :\n");
					poly_print(quo1);
					printf("quo2 :\n");
					poly_print(quo2);
					printf("rem1 :\n");
					poly_print(rem1);
					printf("rem2 :\n");
					poly_print(rem2);
					printf("verif :\n");
					poly_print(verif);
					fprintf(stderr, "1 is not inversible!\n");
					exit(EXIT_FAILURE);
				}
				poly_free_multi(7, pa, pb, quo1, quo2, rem1, rem2, verif);
			}
		}
	}
	return bon;
}

bool comp_gcd(int qqt)
{
	bool bon = true;
	for (int dega = -1; dega < n; dega++)
	{
		for (int degb = -1; degb < n; degb++)
		{
			if (dega != -1 || degb != -1)
			{
				for (int i = 0; i < qqt; i++)
				{
					poly pa = poly_new_rand(dega);
					poly pb = poly_new_rand(degb);
					poly d1 = poly_new();
					poly u1 = poly_new();
					poly v1 = poly_new();
					poly d2 = poly_new();
					poly u2 = poly_new();
					poly v2 = poly_new();
					poly_xgcd(d1, u1, v1, pa, pb);
					//poly_xgcd(d2, u2, v2, pa, pb);
					//poly_fast_xgcd(d1, u1, v1, pa, pb);
					poly_fast_xgcd(d2, u2, v2, pa, pb);
					if (poly_leading_coeff(d1) == 0 || poly_leading_coeff(d2) == 0)
					{
						printf("ca fait n'importe quoi!\n");
						poly_print(pa);
						poly_print(pb);
						poly_print(d1);
						poly_print(u1);
						poly_print(v1);
						poly_print(d2);
						poly_print(u2);
						poly_print(v2);
					}
					int coeff = zp_inv(poly_leading_coeff(d1));
					poly_mul_scalar(d1, coeff, d1);
					poly_mul_scalar(u1, coeff, u1);
					poly_mul_scalar(v1, coeff, v1);
					coeff = zp_inv(poly_leading_coeff(d2));
					poly_mul_scalar(d2, coeff, d2);
					poly_mul_scalar(u2, coeff, u2);
					poly_mul_scalar(v2, coeff, v2);
					poly verif1 = poly_new();
					poly verif2 = poly_new();
					poly prod1 = poly_new();
					poly prod2 = poly_new();
					poly_mul(prod1, u1, pa);
					poly_mul(prod2, v1, pb);
					poly_add(verif1, prod1, prod2);
					poly_mul(prod1, u2, pa);
					poly_mul(prod2, v2, pb);
					poly_add(verif2, prod1, prod2);
					if (!poly_equal(d1, d2) || !poly_equal(d1, verif1) || !poly_equal(d2, verif2))
					{
						bon = false;
						printf("resultats distincts xgcd et fast_xgcd\n");
						poly_print(pa);
						poly_print(pb);
						poly_print(d1);
						poly_print(u1);
						poly_print(v1);
						poly_print(d2);
						poly_print(u2);
						poly_print(v2);
					}
					poly_free_multi(12, pa, pb, d1, u1, v1, d2, u2, v2, verif1, verif2, prod1, prod2);
				}
			}
		}
	}
	return bon;
}