#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <sys\time.h>

#include "field_settings.h"
#include "finite_field.h"
#include "array.h"
#include "polynomial.h"
#include "rs_code.h"

int p = 0;
int q = 0;
int n = 0;
zp_t omega = 0;
array omegas = NULL;
array inverses = NULL;

long long timeInMilliseconds(void) {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000)+(tv.tv_usec/1000);
}

int time_encode(int block_length, int message_length, int qty, array *messages)
{
	long long start = timeInMilliseconds();
	for (int i = 0; i < qty; i++)
		rs_encode_2(block_length, message_length, messages[i]);
	return timeInMilliseconds() - start;
}

int time_fast_encode(int block_length, int message_length, int qty, array *messages)
{
	long long start = timeInMilliseconds();
	for (int i = 0; i < qty; i++)
		rs_fast_encode(block_length, message_length, messages[i]);
	return timeInMilliseconds() - start;
}

int time_decode(poly g_0, int block_length, int message_length, int qty, array *receiveds)
{
	long long start = timeInMilliseconds();
	for (int i = 0; i < qty; i++)
		rs_decode_2(g_0, block_length, message_length, receiveds[i]);
	return timeInMilliseconds() - start;
}

int time_fast_decode(poly g_0, int block_length, int message_length, int qty, array *receiveds)
{
	long long start = timeInMilliseconds();
	for (int i = 0; i < qty; i++)
		rs_fast_decode(g_0, block_length, message_length, receiveds[i]);
	return timeInMilliseconds() - start;
}

void comp_time_encode(int p_, int block_length, int message_length, int qty, int *time1, int *time2)
{
	field_settings_reset(p_);
	array *messages = (array *)malloc(qty * sizeof(array));
	assert(messages);
	for (int i = 0; i < qty; i++)
		messages[i] = array_new_rand(message_length);
	*time1 = time_encode(block_length, message_length, qty, messages);
	*time2 = time_fast_encode(block_length, message_length, qty, messages);
	free(messages);
}

void comp_time_decode(int p_, int block_length, int message_length, int qty, int *time1, int *time2)
{
	field_settings_reset(p_);
	poly g_0 = poly_new();
	rs_g_0_fourier(g_0, block_length);
	array *receiveds = (array *)malloc(qty * sizeof(array));
	assert(receiveds);
	for (int i = 0; i < qty; i++)
		receiveds[i] = rs_encode_2(block_length, message_length, array_new_rand(message_length));
	*time1 = time_decode(g_0, block_length, message_length, qty, receiveds);
	*time2 = time_fast_decode(g_0, block_length, message_length, qty, receiveds);
	free(receiveds);
	poly_free(g_0);
}

void write_results(int ***results, int nb_prime, int *block_lengths, int max_block_length, int primes[])
{
	FILE *fptr_encode;
	FILE *fptr_decode;
	fptr_encode = fopen("results_encode.csv", "w");
	assert(fptr_encode);
	fptr_decode = fopen("results_decode.csv", "w");
	assert(fptr_decode);
	fprintf(fptr_encode, "\"message length\"");
	fprintf(fptr_decode, "\"message length\"");
	for (int i = 0; i < nb_prime; i++)
	{
		fprintf(fptr_encode, ",\"encode p = %d\",\"fast encode p = %d\"", primes[i], primes[i]);
		fprintf(fptr_decode, ",\"decode p = %d\",\"fast decode p = %d\"", primes[i], primes[i]);
	}
	fprintf(fptr_encode, "\n");
	fprintf(fptr_decode, "\n");
	for (int message_length = 1; message_length < max_block_length; message_length++)
	{
		fprintf(fptr_encode, "%d", message_length);
		fprintf(fptr_decode, "%d", message_length);
		for (int i = 0; i < nb_prime; i++)
		{
			if (block_lengths[i] > message_length)
			{
				fprintf(fptr_encode, ",%d,%d", results[i][message_length][0], results[i][message_length][1]);
				fprintf(fptr_decode, ",%d,%d", results[i][message_length][2], results[i][message_length][3]);
			}
			else
			{
				fprintf(fptr_encode, ",,");
				fprintf(fptr_decode, ",,");
			}
		}
		fprintf(fptr_encode, "\n");
		fprintf(fptr_decode, "\n");
	}
	fclose(fptr_encode);
	fclose(fptr_decode);
}

int main(void)
{
	field_settings_set(193);
	int primes[32] = {3, 5, 7, 11, 13, 17, 19, 29, 37, 41, 73, 97, 113, 193, 257, 449, 577, 641, 769, 1153, 12289, 18433, 40961, 65537, 114689, 147457, 163841, 786433, 1179649, 7340033, 167772161, 469762049};
	int nb_prime = 20;
	int *block_lengths = (int *)malloc(nb_prime * sizeof(int));
	assert(block_lengths);
	int max_block_length = 0;
	for (int i = 0; i < nb_prime; i++)
	{
		block_lengths[i] = 1;
		while (((primes[i] - 1) % block_lengths[i]) == 0)
			block_lengths[i] *= 2;
		block_lengths[i] /= 2;
		if (block_lengths[i] > max_block_length)
			max_block_length = block_lengths[i];
	}
	int ***results = (int ***)malloc(nb_prime * sizeof(int **));
	assert(results);
	int time1 = 0;
	int time2 = 0;
	for (int i = 0; i < nb_prime; i++)
	{
		printf("prime %d / %d, p = %d, n = %d\n", (i + 1), nb_prime, primes[i], block_lengths[i]);
		results[i] = (int **)malloc(block_lengths[i] * sizeof(int *));
		assert(results[i]);
		for (int message_length = 1; message_length < block_lengths[i]; message_length++)
		{
			printf("message_length = %d / %d\n", message_length, block_lengths[i] - 1);
			results[i][message_length] = (int *)malloc(4 * sizeof(int));
			assert(results[i][message_length]);
			comp_time_encode(primes[i], block_lengths[i], message_length, 20, &time1, &time2);
			results[i][message_length][0] = time1;
			results[i][message_length][1] = time2;
			comp_time_decode(primes[i], block_lengths[i], message_length, 20, &time1, &time2);
			results[i][message_length][2] = time1;
			results[i][message_length][3] = time2;
			printf("\033[F\033[K");
		}
	}
	write_results(results, nb_prime, block_lengths, max_block_length, primes);
	for (int i = 0; i < nb_prime; i++)
	{
		for (int message_length = 1; message_length < block_lengths[i]; message_length++)
			free(results[i][message_length]);
		free(results[i]);
	}
	free(block_lengths);
	free(results);
	field_settings_free();
	return EXIT_SUCCESS;
}