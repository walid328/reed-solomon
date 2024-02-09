#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>

#include "../src/field_settings.h"
#include "../src/finite_field.h"
#include "../src/zp_array.h"
#include "../src/polynomial.h"
#include "../src/rs_code.h"

int p = 0;
int q = 0;
int n = 0;
zp_t omega = 0;
zp_array omegas = NULL;
zp_array inverses = NULL;

int rand_range(int min, int max)
{
	int ret = (rand() << 24) ^ (rand() << 16) ^ (rand() << 8) ^ rand();
	ret %= (max - min + 1);
	if (ret < 0)
		ret += max - min + 1;
	return ret + min;
}

int maximum(int a, int b)
{
	return a > b ? a : b;
}

long long timeInMilliseconds(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (((long long)tv.tv_sec) * 1000) + (tv.tv_usec / 1000);
}

int time_mul(int qty, poly *op_1s, poly *op_2s)
{
	long long start = timeInMilliseconds();
	poly rop = poly_new();
	for (int i = 0; i < qty; i++)
		poly_mul(rop, op_1s[i], op_2s[i]);
	poly_free(rop);
	return timeInMilliseconds() - start;
}

int time_fast_mul(int qty, poly *op_1s, poly *op_2s)
{
	long long start = timeInMilliseconds();
	poly rop = poly_new();
	for (int i = 0; i < qty; i++)
		poly_fast_mul(rop, op_1s[i], op_2s[i]);
	poly_free(rop);
	return timeInMilliseconds() - start;
}

int time_euc_div(int qty, poly *op_1s, poly *op_2s)
{
	long long start = timeInMilliseconds();
	poly quo = poly_new();
	poly rem = poly_new();
	for (int i = 0; i < qty; i++)
		poly_euc_div(quo, rem, op_1s[i], op_2s[i]);
	poly_free_multi(2, quo, rem);
	return timeInMilliseconds() - start;
}

int time_fast_euc_div(int qty, poly *op_1s, poly *op_2s)
{
	long long start = timeInMilliseconds();
	poly quo = poly_new();
	poly rem = poly_new();
	for (int i = 0; i < qty; i++)
		poly_fast_euc_div(quo, rem, op_1s[i], op_2s[i]);
	poly_free_multi(2, quo, rem);
	return timeInMilliseconds() - start;
}

int time_xgcd(int qty, poly *op_1s, poly *op_2s)
{
	long long start = timeInMilliseconds();
	poly d = poly_new();
	poly u = poly_new();
	poly v = poly_new();
	for (int i = 0; i < qty; i++)
		poly_xgcd(d, u, v, op_1s[i], op_2s[i]);
	poly_free_multi(3, d, u, v);
	return timeInMilliseconds() - start;
}

int time_fast_xgcd(int qty, poly *op_1s, poly *op_2s)
{
	long long start = timeInMilliseconds();
	poly d = poly_new();
	poly u = poly_new();
	poly v = poly_new();
	for (int i = 0; i < qty; i++)
		poly_fast_xgcd(d, u, v, op_1s[i], op_2s[i]);
	poly_free_multi(3, d, u, v);
	return timeInMilliseconds() - start;
}

int time_encode(int block_length, int *message_lengths, int qty, zp_array *messages)
{
	long long start = timeInMilliseconds();
	for (int i = 0; i < qty; i++)
		zp_array_free(rs_encode_2(block_length, message_lengths[i], messages[i]));
	return timeInMilliseconds() - start;
}

int time_fast_encode(int block_length, int *message_lengths, int qty, zp_array *messages)
{
	long long start = timeInMilliseconds();
	for (int i = 0; i < qty; i++)
		zp_array_free(rs_fast_encode(block_length, message_lengths[i], messages[i]));
	return timeInMilliseconds() - start;
}

int time_decode(poly g_0, int block_length, int *message_lengths, int qty, zp_array *receiveds)
{
	long long start = timeInMilliseconds();
	for (int i = 0; i < qty; i++)
		zp_array_free(rs_decode_2(g_0, block_length, message_lengths[i], receiveds[i]));
	return timeInMilliseconds() - start;
}

int time_fast_decode(poly g_0, int block_length, int *message_lengths, int qty, zp_array *receiveds)
{
	long long start = timeInMilliseconds();
	for (int i = 0; i < qty; i++)
		zp_array_free(rs_fast_decode(g_0, block_length, message_lengths[i], receiveds[i]));
	return timeInMilliseconds() - start;
}

void comp_time_mul(int qty, int *time1, int *time2)
{
	poly *op_1s = (poly *)malloc(qty * sizeof(poly));
	assert(op_1s);
	poly *op_2s = (poly *)malloc(qty * sizeof(poly));
	assert(op_2s);
	for (int i = 0; i < qty; i++)
	{
		op_1s[i] = poly_new_rand(rand_range(-1, 2 * n));
		op_2s[i] = poly_new_rand(rand_range(-1, 2 * n - poly_deg(op_1s[i]) - 1));
	}
	*time1 = time_mul(qty, op_1s, op_2s);
	*time2 = time_fast_mul(qty, op_1s, op_2s);
	for (int i = 0; i < qty; i++)
		poly_free_multi(2, op_1s[i], op_2s[i]);
	free(op_1s);
	free(op_2s);
}

void comp_time_euc_div(int qty, int *time1, int *time2)
{
	poly *op_1s = (poly *)malloc(qty * sizeof(poly));
	assert(op_1s);
	poly *op_2s = (poly *)malloc(qty * sizeof(poly));
	assert(op_2s);
	for (int i = 0; i < qty; i++)
	{
		op_1s[i] = poly_new_rand(rand_range(-1, 2 * n - 1));
		op_2s[i] = poly_new_rand(rand_range(maximum(0, poly_deg(op_1s[i]) + 1 - n), n));
	}
	*time1 = time_euc_div(qty, op_1s, op_2s);
	*time2 = time_fast_euc_div(qty, op_1s, op_2s);
	for (int i = 0; i < qty; i++)
		poly_free_multi(2, op_1s[i], op_2s[i]);
	free(op_1s);
	free(op_2s);
}

void comp_time_xgcd(int qty, int *time1, int *time2)
{
	poly *op_1s = (poly *)malloc(qty * sizeof(poly));
	assert(op_1s);
	poly *op_2s = (poly *)malloc(qty * sizeof(poly));
	assert(op_2s);
	for (int i = 0; i < qty; i++)
	{
		op_1s[i] = poly_new_rand(rand_range(-1, n));
		op_2s[i] = poly_new_rand(rand_range(poly_deg(op_1s[i]) == -1 ? 0 : -1, n));
	}
	*time1 = time_xgcd(qty, op_1s, op_2s);
	*time2 = time_fast_xgcd(qty, op_1s, op_2s);
	for (int i = 0; i < qty; i++)
		poly_free_multi(2, op_1s[i], op_2s[i]);
	free(op_1s);
	free(op_2s);
}

void comp_time_encode(int qty, int *time1, int *time2)
{
	zp_array *messages = (zp_array *)malloc(qty * sizeof(zp_array));
	assert(messages);
	int *message_lengths = (int *)malloc(qty * sizeof(int));
	assert(message_lengths);
	for (int i = 0; i < qty; i++)
	{
		message_lengths[i] = rand_range(1, n - 1);
		messages[i] = zp_array_new_rand(message_lengths[i]);
	}
	*time1 = time_encode(n, message_lengths, qty, messages);
	*time2 = time_fast_encode(n, message_lengths, qty, messages);
	for (int i = 0; i < qty; i++)
		zp_array_free(messages[i]);
	free(messages);
	free(message_lengths);
}

void comp_time_decode(int qty, int *time1, int *time2)
{
	poly g_0 = poly_new();
	rs_g_0_fourier(g_0, n);
	zp_array *receiveds = (zp_array *)malloc(qty * sizeof(zp_array));
	assert(receiveds);
	int *message_lengths = (int *)malloc(qty * sizeof(int));
	assert(message_lengths);
	for (int i = 0; i < qty; i++)
	{
		message_lengths[i] = rand_range(1, n - 1);
		zp_array message = zp_array_new_rand(message_lengths[i]);
		receiveds[i] = rs_encode_2(n, message_lengths[i], message);
		zp_array_add_errors(receiveds[i], n, rand_range(0, (n - message_lengths[i]) / 2));
		zp_array_free(message);
	}
	*time1 = time_decode(g_0, n, message_lengths, qty, receiveds);
	*time2 = time_fast_decode(g_0, n, message_lengths, qty, receiveds);
	for (int i = 0; i < qty; i++)
		zp_array_free(receiveds[i]);
	free(receiveds);
	free(message_lengths);
	poly_free(g_0);
}

void write_results(int **results, int nb_prime, int primes[])
{
	FILE *fptr_results;
	fptr_results = fopen("results_perf/results.csv", "w");
	assert(fptr_results);
	fprintf(fptr_results, "\"prime p\"");
	fprintf(fptr_results, ",\"mul\",\"fast_mul\"");
	fprintf(fptr_results, ",\"euc_div\",\"fast_euc_div\"");
	fprintf(fptr_results, ",\"xgcd\",\"fast_xgcd\"");
	fprintf(fptr_results, ",\"encode\",\"fast_encode\"");
	fprintf(fptr_results, ",\"decode\",\"fast_decode\"");
	fprintf(fptr_results, "\n");
	for (int i = 0; i < nb_prime; i++)
	{
		fprintf(fptr_results, "%d", primes[i]);
		for (int j = 0; j < 10; j++)
			fprintf(fptr_results, ",%d", results[i][j]);
		fprintf(fptr_results, "\n");
	}
	fclose(fptr_results);
}

int main(void)
{
	srand(time(NULL));
	field_settings_set(193);
	// Here primes are sorted by n increasing, n is the biggest power of 2 dividing p-1.
	int primes[32] = {3, 7, 11, 19, 5, 13, 29, 37, 41, 73, 17, 113, 97, 193, 449, 577, 641, 1153, 257, 769, 18433, 12289, 40961, 114689, 147457, 163841, 65537, 1179649, 786433, 7340033, 167772161, 469762049};
	// Here primes are sorted in ascending order.
	// int primes[32] = {3, 5, 7, 11, 13, 17, 19, 29, 37, 41, 73, 97, 113, 193, 257, 449, 577, 641, 769, 1153, 12289, 18433, 40961, 65537, 114689, 147457, 163841, 786433, 1179649, 7340033, 167772161, 469762049};
	int nb_prime = 23;
	int **results = (int **)malloc(nb_prime * sizeof(int *));
	assert(results);
	int time1 = 0;
	int time2 = 0;
	for (int i = 0; i < nb_prime; i++)
	{
		field_settings_set(primes[i]);
		printf("prime %d / %d, p = %d, n = %d\n", (i + 1), nb_prime, primes[i], n);
		results[i] = (int *)malloc(10 * sizeof(int));
		assert(results[i]);
		int qty = 100;
		comp_time_mul(qty, &time1, &time2);
		results[i][0] = time1;
		results[i][1] = time2;
		comp_time_euc_div(qty, &time1, &time2);
		results[i][2] = time1;
		results[i][3] = time2;
		comp_time_xgcd(qty, &time1, &time2);
		results[i][4] = time1;
		results[i][5] = time2;
		comp_time_encode(qty, &time1, &time2);
		results[i][6] = time1;
		results[i][7] = time2;
		comp_time_decode(qty, &time1, &time2);
		results[i][8] = time1;
		results[i][9] = time2;
	}
	write_results(results, nb_prime, primes);
	for (int i = 0; i < nb_prime; i++)
		free(results[i]);
	free(results);
	field_settings_free();
	return EXIT_SUCCESS;
}