#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "finite_field.h"

struct polynomial
{
	int degree;
	int *coefficients;
};

int p = 7;

int size_int_str(int n)
{
	int ret = 1;
	if (n < 0)
	{
		ret = 2;
		n = -n;
	}
	for (; n >= 10; n /= 10)
		ret++;
	return ret;
}

void str_add_int(char *string, int *index, int n)
{
	if (n == 0)
	{
		string[*index] = '0';
		*index += 1;
		return;
	}
	if (n < 0)
	{
		string[*index] = '-';
		*index += 1;
		n = -n;
	}
	int size_n = size_int_str(n);
	for (int i = size_n - 1; i >= 0; i--, n /= 10)
		string[*index + i] = n % 10;
	*index += size_n;
}

char *str_poly_iso_length(poly *q)
{
	if (q->degree == -1)
	{
		char *poly_string = (char *)malloc(2 * sizeof(char));
		assert(poly_string != NULL);
		poly_string[0] = '0';
		poly_string[1] = '\0';
		return poly_string;
	}
	int size_p = size_int_str(p);
	int size_d = size_int_str(q->degree);
	char *poly_string = (char *)malloc(((size_p + size_d + 6) * (q->degree + 1) + 1) * sizeof(char));
	assert(poly_string != NULL);
	int index = 0;
	for (int i = q->degree; i >= 0; i--)
	{
		if (i != q->degree)
		{
			poly_string[index + 0] = ' ';
			poly_string[index + 1] = '+';
			poly_string[index + 2] = ' ';
			index += 3;
		}
		int size_coeff = size_int_str(q->coefficients[i]);
		for (int j = 0; j < size_p - size_coeff; j++)
			poly_string[index + j] = ' ';
		index += size_p - size_coeff;
		str_add_int(poly_string, &index, q->coefficients[i]);
		poly_string[index + 0] = '*';
		poly_string[index + 1] = 'X';
		poly_string[index + 2] = '^';
		index += 3;
		str_add_int(poly_string, &index, i);
		int size_exponant = size_int_str(i);
		for (int j = 0; j < size_d - size_exponant; j++)
			poly_string[index + j] = ' ';
		index += size_d - size_exponant;
	}
	poly_string[index] = '\0';
	return poly_string;
}

char *str_poly_min(poly *q)
{
	if (q->degree == -1)
	{
		char *poly_string = (char *)malloc(2 * sizeof(char));
		assert(poly_string != NULL);
		poly_string[0] = '0';
		poly_string[1] = '\0';
		return poly_string;
	}
	int size_p = size_int_str(p);
	int size_d = size_int_str(q->degree);
	char *poly_string = (char *)malloc(((size_p + size_d + 6) * (q->degree + 1) + 1) * sizeof(char));
	assert(poly_string != NULL);
	int index = 0;
	for (int i = q->degree; i >= 0; i--)
	{
		if (q->coefficients[i] != 0)
		{
			if (i != q->degree)
			{
				poly_string[index + 0] = ' ';
				poly_string[index + 1] = '+';
				poly_string[index + 2] = ' ';
				index += 3;
			}
			str_add_int(poly_string, &index, q->coefficients[i]);
			if (i > 0)
			{
				poly_string[index + 0] = '*';
				poly_string[index + 1] = 'X';
				index += 2;
				if (i > 1)
				{
					poly_string[index + 0] = '^';
					index += 1;
					str_add_int(poly_string, &index, i);
				}
			}
		}
	}
	poly_string[index] = '\0';
	return poly_string;
}

void print_poly_iso_length(poly *q)
{
	char *poly_string = str_poly_iso_length(q);
	printf("%s\n", poly_string);
	free(poly_string);
}

void print_poly_min(poly *q)
{
	char *poly_string = str_poly_min(q);
	printf("%s\n", poly_string);
	free(poly_string);
}

void print_poly(poly *q)
{
	if (q->degree == -1)
		printf("0");
	for (int i = q->degree; i >= 0; i--)
	{
		if (i != q->degree)
			printf(" + ");
		printf("%d", q->coefficients[i]);
		if (i == 1)
			printf("*X");
		if (i > 1)
			printf("*X^%d", i);
	}
	printf("\n");
}

poly *new_poly(void)
{
	poly *q = (poly *)malloc(sizeof(poly));
	assert(q);
	q->degree = -1;
	q->coefficients = NULL;
	return q;
}

void set_poly(poly *q, int degree, int *coefficients)
{
	q->degree = degree;
	q->coefficients = coefficients;
}

void free_poly(poly *q)
{
	free(q);
}

poly *add_poly(poly *p1, poly *p2)
{
	poly *sum = new_poly();
	if (p1->degree > p2->degree)
	{
		int *coefficients = (int *)malloc(p1->degree * sizeof(int));
		assert(coefficients);
		int i = 0;
		for (; i <= p2->degree; i++)
			coefficients[i] = (p1->coefficients[i] + p2->coefficients[i]) % p;
		for (; i <= p1->degree; i++)
			coefficients[i] = p1->coefficients[i];
		set_poly(sum, p1->degree, coefficients);
	}
	else if (p1->degree < p2->degree)
	{
		int *coefficients = (int *)malloc(p1->degree * sizeof(int));
		assert(coefficients);
		int i = 0;
		for (; i <= p1->degree; i++)
			coefficients[i] = (p1->coefficients[i] + p2->coefficients[i]) % p;
		for (; i <= p2->degree; i++)
			coefficients[i] = p2->coefficients[i];
		set_poly(sum, p2->degree, coefficients);
	}
	else if (p1->degree == -1)
	{
		int *coefficient = (int *)calloc(1, sizeof(int));
		assert(coefficient);
		set_poly(sum, -1, coefficient);
	}
	else
	{
		int i;
		for (i = p1->degree; i >= 0 && (p1->coefficients[i] + p2->coefficients[i]) % p == 0; i--)
			;
		if (i == -1)
		{
			int *coefficient = (int *)calloc(1, sizeof(int));
			assert(coefficient);
			set_poly(sum, -1, coefficient);
		}
		else
		{
			int degree = i;
			int *coefficients = (int *)malloc(i * sizeof(int));
			assert(coefficients);
			for (; i >= 0; i--)
				coefficients[i] = (p1->coefficients[i] + p2->coefficients[i]) % p;
			set_poly(sum, degree, coefficients);
		}
	}
	return sum;
}

poly *mul_poly(poly *p1, poly *p2)
{
	poly *prod = new_poly();
	if (p1->degree == -1 || p2->degree == 0)
	{
		int *coefficient = (int *)calloc(1, sizeof(int));
		set_poly(prod, -1, coefficient);
	}
	else
	{
		int degree = p1->degree + p2->degree;
		int *coefficients = (int *)calloc(degree, sizeof(int));
		assert(coefficients);
		for (int i = 0; i <= p1->degree; i++)
			for (int j = 0; j <= p2->degree; j++)
				coefficients[i + j] = (coefficients[i + j] + (p1->coefficients[i] * p2->coefficients[j])) % p;
		set_poly(prod, degree, coefficients);
	}
	return prod;
}