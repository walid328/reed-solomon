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

// Tests to ensure that fast operations works the same.

int p = 0;
int q = 0;
int n = 0;
zp_t omega = 0;
array omegas = NULL;
array inverses = NULL;

void usage(int argc, char *argv[])
{
	fprintf(stderr, "\e[1;31mFailure: \e[0mNot enough parameters!\n");
	fprintf(stderr, "\e[1mUsage:\e[0m %s <p> <qty>\n", argv[0]);
	exit(EXIT_FAILURE);
}

void print_result(bool res, char *test_name)
{
	if (res)
		printf("Test \"%s\" finished: \e[1;32mSUCCESS\e[0m\n", test_name);
	else
		printf("Test \"%s\" finished: \e[1;31mFAILURE\e[0m\n", test_name);
}

bool comp_dft(int qty)
{
	for (int dega = -1; dega < n; dega++)
	{
		for (int i = 0; i < qty; i++)
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
					printf("Execution stopped for values:\n");
					poly_print(pa);
					array_print(eval1, d);
					array_print(eval2, d);
					return false;
				}
				array_free(eval1);
				array_free(eval2);
				d *= 2;
			}
			poly_free(pa);
		}
	}
	return true;
}

bool comp_inv_dft(int qty)
{
	for (int d = 1; d <= n; d *= 2)
	{
		for (int i = 0; i < qty; i++)
		{
			array eval = array_new_rand(d);
			poly pa = poly_new();
			poly pb = poly_new();
			poly_inv_dft(pa, eval, d);
			poly_inv_fft(pb, eval, d);
			if (!poly_equal(pa, pb))
			{
				printf("Execution stopped for values:\n");
				array_print(eval, d);
				poly_print(pa);
				poly_print(pb);
				return false;
			}
			array_free(eval);
			poly_free_multi(2, pa, pb);
		}
	}
	return true;
}

bool comp_mul(int qty)
{
	for (int dega = -1; dega <= 2 * n; dega++)
	{
		for (int degb = -1; dega + degb + 1 <= 2 * n; degb++)
		{
			for (int i = 0; i < qty; i++)
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
					printf("Execution stopped for values:\n");
					printf("%d\n", poly_equal(prod1, prod2));
					printf("%d\n", poly_equal(prod1, prod3));
					printf("%d\n", poly_equal(prod1, prod4));
					poly_print(pa);
					poly_print(pb);
					poly_print(prod1);
					poly_print(prod2);
					poly_print(prod3);
					poly_print(prod4);
					return false;
				}
				poly_free_multi(6, pa, pb, prod1, prod2, prod3, prod4);
			}
		}
	}
	return true;
}

bool comp_euc_div(int qty)
{
	for (int dega = -1; dega <= 2 * n - 1; dega++)
	{
		int min_degb = 0;
		if (dega + 1 - n > min_degb)
			min_degb = dega + 1 - n;
		for (int degb = min_degb; degb <= n; degb++)
		{
			for (int i = 0; i < qty; i++)
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
					printf("Execution stopped for values:\n");
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
					return false;
				}
				poly_free_multi(7, pa, pb, quo1, quo2, rem1, rem2, verif);
			}
		}
	}
	return true;
}

bool comp_xgcd(int qty)
{
	for (int dega = -1; dega <= n; dega++)
	{
		for (int degb = -1; degb <= n; degb++)
		{
			if (dega != -1 || degb != -1)
			{
				for (int i = 0; i < qty; i++)
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
					poly_fast_xgcd(d2, u2, v2, pa, pb);
					zp_t coeff = zp_inv(poly_leading_coeff(d1));
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
						printf("Execution stopped for values:\n");
						poly_print(pa);
						poly_print(pb);
						poly_print(d1);
						poly_print(u1);
						poly_print(v1);
						poly_print(d2);
						poly_print(u2);
						poly_print(v2);
						return false;
					}
					poly_free_multi(12, pa, pb, d1, u1, v1, d2, u2, v2, verif1, verif2, prod1, prod2);
				}
			}
		}
	}
	return true;
}

bool comp_xgcd_partial(int qty)
{
	for (int dega = -1; dega <= n; dega++)
	{
		for (int degb = -1; degb <= n; degb++)
		{
			// printf("dega = %d , degb = %d\n", dega, degb);
			if (dega != -1 || degb != -1)
			{
				for (int limit = 0; limit <= dega + 1 && limit <= degb + 1; limit++)
				{
					for (int i = 0; i < qty; i++)
					{
						poly pa = poly_new_rand(dega);
						poly pb = poly_new_rand(degb);
						poly d1 = poly_new();
						poly u1 = poly_new();
						poly v1 = poly_new();
						poly d2 = poly_new();
						poly u2 = poly_new();
						poly v2 = poly_new();
						poly_xgcd_partial(d1, u1, v1, pa, pb, limit);
						poly_fast_xgcd_partial(d2, u2, v2, pa, pb, limit);
						zp_t coeff = zp_inv(poly_leading_coeff(d1));
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
							printf("Execution stopped for values:\n");
							printf("limit = %d\n", limit);
							poly_print(pa);
							poly_print(pb);
							poly_print(d1);
							poly_print(u1);
							poly_print(v1);
							poly_print(d2);
							poly_print(u2);
							poly_print(v2);
							return false;
						}
						poly_free_multi(12, pa, pb, d1, u1, v1, d2, u2, v2, verif1, verif2, prod1, prod2);
					}
				}
			}
		}
	}
	return true;
}

bool comp_encode(int qty)
{
	for (int block_length = 2; block_length <= n; block_length *= 2)
	{
		for (int message_length = 1; message_length < block_length; message_length++)
		{
			for (int i = 0; i < qty; i++)
			{
				array message = array_new_rand(message_length);
				array codeword1 = rs_encode_2(block_length, message_length, message);
				array codeword2 = rs_fast_encode(block_length, message_length, message);
				if (!array_equal(codeword1, codeword2, block_length))
				{
					printf("Execution stopped for values:\n");
					array_print(message, message_length);
					array_print(codeword1, block_length);
					array_print(codeword2, block_length);
					return false;
				}
				array_free(message);
				array_free(codeword1);
				array_free(codeword2);
			}
		}
	}
	return true;
}

bool comp_decode(int qty)
{
	for (int block_length = 2; block_length <= n; block_length *= 2)
	{
		poly g_0 = poly_new();
		rs_g_0_fourier(g_0, block_length);
		for (int message_length = 1; message_length < block_length; message_length++)
		{
			for (int number_errors = 0; number_errors <= (block_length - message_length) / 2; number_errors++)
			{
				for (int i = 0; i < qty; i++)
				{
					array message = array_new_rand(message_length);
					array received = rs_fast_encode(block_length, message_length, message);
					array_add_errors(received, block_length, number_errors);
					array decode1 = rs_decode_2(g_0, block_length, message_length, received);
					array decode2 = rs_fast_decode(g_0, block_length, message_length, received);
					if (!array_equal(decode1, message, message_length) || !array_equal(decode1, decode2, message_length))
					{
						printf("Execution stopped for values:\n");
						array_print(message, message_length);
						array_print(received, block_length);
						array_print(decode1, message_length);
						array_print(decode2, message_length);
						return false;
					}
					array_free(message);
					array_free(received);
					array_free(decode1);
					array_free(decode2);
				}
			}
		}
		poly_free(g_0);
	}
	return true;
}

int main(int argc, char *argv[])
{
	if (argc < 3)
		usage(argc, argv);

	int p_ = atoi(argv[1]);
	field_settings_set(p_);

	int qty = 1;

	bool ok;
	if (strcmp("dft", argv[2]) == 0)
		ok = comp_dft(qty);
	else if (strcmp("inv_dft", argv[2]) == 0)
		ok = comp_inv_dft(qty);
	else if (strcmp("mul", argv[2]) == 0)
		ok = comp_mul(qty);
	else if (strcmp("euc_div", argv[2]) == 0)
		ok = comp_euc_div(qty);
	else if (strcmp("xgcd", argv[2]) == 0)
		ok = comp_xgcd(qty);
	else if (strcmp("xgcd_partial", argv[2]) == 0)
		ok = comp_xgcd_partial(qty);
	else if (strcmp("encode", argv[2]) == 0)
		ok = comp_encode(qty);
	else if (strcmp("decode", argv[2]) == 0)
		ok = comp_decode(qty);
	else if (strcmp("compare", argv[2]) == 0 || strcmp("all", argv[2]) == 0 || strcmp("*", argv[2]) == 0)
	{
		bool sub_ok = comp_dft(qty);
		print_result(sub_ok, "comp_dft");
		ok = sub_ok;
		sub_ok = comp_inv_dft(qty);
		print_result(sub_ok, "comp_inv_dft");
		ok &= sub_ok;
		sub_ok = comp_mul(qty);
		print_result(sub_ok, "comp_mul");
		ok &= sub_ok;
		sub_ok = comp_euc_div(qty);
		print_result(sub_ok, "comp_euc_div");
		ok &= sub_ok;
		sub_ok = comp_xgcd(qty);
		print_result(sub_ok, "comp_xgcd");
		ok &= sub_ok;
		sub_ok = comp_xgcd_partial(qty);
		print_result(sub_ok, "comp_partial");
		ok &= sub_ok;
		sub_ok = comp_encode(qty);
		print_result(sub_ok, "comp_encode");
		ok &= sub_ok;
		sub_ok = comp_decode(qty);
		print_result(sub_ok, "comp_decode");
		ok &= sub_ok;
	}

	else
	{
		fprintf(stderr, "Error: test \"%s\" not found!\n", argv[2]);
		exit(EXIT_FAILURE);
	}

	printf("==> ");
	print_result(ok, argv[2]);

	field_settings_free();
	return EXIT_SUCCESS;
}
