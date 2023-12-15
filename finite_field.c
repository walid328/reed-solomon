#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "finite_field.h"

struct polynomial
{
	int degree;
	int *coeffs;
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
		string[*index + i] = '0' + (n % 10);
	*index += size_n;
}

char *str_poly_iso_length(poly *q)
{
	if (q->degree == -1)
	{
		char *poly_string = (char *)malloc(2 * sizeof(char));
		assert(poly_string);
		poly_string[0] = '0';
		poly_string[1] = '\0';
		return poly_string;
	}
	int size_p = size_int_str(p);
	int size_d = size_int_str(q->degree);
	char *poly_string = (char *)malloc(((size_p + size_d + 6) * (q->degree + 1) + 1) * sizeof(char));
	assert(poly_string);
	int index = 0;
	for (int i = 0; i <= q->degree; i++)
	{
		if (index != 0)
		{
			poly_string[index + 0] = ' ';
			poly_string[index + 1] = '+';
			poly_string[index + 2] = ' ';
			index += 3;
		}
		int size_coeff = size_int_str(q->coeffs[i]);
		for (int j = 0; j < size_p - size_coeff; j++)
			poly_string[index + j] = ' ';
		index += size_p - size_coeff;
		str_add_int(poly_string, &index, q->coeffs[i]);
		if (i == 1)
		{
			poly_string[index + 0] = '*';
			poly_string[index + 1] = 'x';
			index += 2;
		}
		if (i > 1)
		{
			poly_string[index + 0] = '*';
			poly_string[index + 1] = 'x';
			poly_string[index + 2] = '^';
			index += 3;
			str_add_int(poly_string, &index, i);
		}
	}
	poly_string[index] = '\0';
	return poly_string;
}

char *str_poly_min(poly *q)
{
	if (q->degree == -1)
	{
		char *poly_string = (char *)malloc(2 * sizeof(char));
		assert(poly_string);
		poly_string[0] = '0';
		poly_string[1] = '\0';
		return poly_string;
	}
	int size_p = size_int_str(p);
	int size_d = size_int_str(q->degree);
	char *poly_string = (char *)malloc(((size_p + size_d + 6) * (q->degree + 1) + 1) * sizeof(char));
	assert(poly_string);
	int index = 0;
	for (int i = 0; i <= q->degree; i++)
	{
		if (q->coeffs[i] != 0)
		{
			if (index != 0)
			{
				poly_string[index + 0] = ' ';
				poly_string[index + 1] = '+';
				poly_string[index + 2] = ' ';
				index += 3;
			}
			str_add_int(poly_string, &index, q->coeffs[i]);
			if (i > 0)
			{
				poly_string[index + 0] = '*';
				poly_string[index + 1] = 'x';
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
	for (int i = 0; i <= q->degree; i++)
	{
		if (i != 0)
			printf(" + ");
		printf("%d", q->coeffs[i]);
		if (i == 1)
			printf("*x");
		if (i > 1)
			printf("*x^%d", i);
	}
	printf("\n");
}

poly *new_poly(void)
{
	poly *q = (poly *)malloc(sizeof(poly));
	assert(q);
	q->degree = -1;
	q->coeffs = NULL;
	return q;
}

poly *new_poly_0(void)
{
	poly *q = (poly *)malloc(sizeof(poly));
	assert(q);
	q->degree = -1;
	q->coeffs = (int *)calloc(1, sizeof(int));
	assert(q->coeffs);
	return q;
}

poly *new_poly_1(void)
{
	poly *q = (poly *)malloc(sizeof(poly));
	assert(q);
	q->degree = 0;
	q->coeffs = (int *)malloc(1 * sizeof(int));
	assert(q->coeffs);
	q->coeffs[0] = 1;
	return q;
}

poly *new_poly_from_copy(poly *p)
{
	if (p->degree == -1)
		return new_poly_0();
	poly *q = (poly *)malloc(sizeof(poly));
	assert(q);
	q->degree = p->degree;
	q->coeffs = (int *)malloc((q->degree + 1) * sizeof(int));
	assert(q->coeffs);
	for (int i = 0; i <= q->degree; i++)
		q->coeffs[i] = p->coeffs[i];
	return q;
}

void set_poly(poly *q, int degree, int *coeffs)
{
	q->degree = degree;
	q->coeffs = coeffs;
}

void free_poly(poly *q)
{
	assert(q);
	free(q);
}

void free_poly_full(poly *q)
{
	assert(q);
	assert(q->coeffs);
	free(q->coeffs);
	free(q);
}

poly *add_poly(poly *p1, poly *p2)
{
	poly *sum = new_poly();
	if (p1->degree > p2->degree)
	{
		int *coeffs = (int *)malloc((p1->degree + 1) * sizeof(int));
		assert(coeffs);
		int i = 0;
		for (; i <= p2->degree; i++)
			coeffs[i] = ((p1->coeffs[i] + p2->coeffs[i]) < p) ? (p1->coeffs[i] + p2->coeffs[i]) : (p1->coeffs[i] + p2->coeffs[i] - p);
		for (; i <= p1->degree; i++)
			coeffs[i] = p1->coeffs[i];
		set_poly(sum, p1->degree, coeffs);
	}
	else if (p1->degree < p2->degree)
	{
		int *coeffs = (int *)malloc((p2->degree + 1) * sizeof(int));
		assert(coeffs);
		int i = 0;
		for (; i <= p1->degree; i++)
			coeffs[i] = ((p1->coeffs[i] + p2->coeffs[i]) < p) ? (p1->coeffs[i] + p2->coeffs[i]) : (p1->coeffs[i] + p2->coeffs[i] - p);
		for (; i <= p2->degree; i++)
			coeffs[i] = p2->coeffs[i];
		set_poly(sum, p2->degree, coeffs);
	}
	else if (p1->degree == -1)
	{
		int *coeff = (int *)calloc(1, sizeof(int));
		assert(coeff);
		set_poly(sum, -1, coeff);
	}
	else
	{
		int i;
		for (i = p1->degree; i >= 0 && (p1->coeffs[i] + p2->coeffs[i]) % p == 0; i--)
			;
		if (i == -1)
		{
			int *coeff = (int *)calloc(1, sizeof(int));
			assert(coeff);
			set_poly(sum, -1, coeff);
		}
		else
		{
			int degree = i;
			int *coeffs = (int *)malloc((degree + 1) * sizeof(int));
			assert(coeffs);
			for (; i >= 0; i--)
				coeffs[i] = ((p1->coeffs[i] + p2->coeffs[i]) < p) ? (p1->coeffs[i] + p2->coeffs[i]) : (p1->coeffs[i] + p2->coeffs[i] - p);
			set_poly(sum, degree, coeffs);
		}
	}
	return sum;
}

poly *substract_poly(poly *p1, poly *p2)
{
	poly *dif = new_poly();
	if (p1->degree > p2->degree)
	{
		int *coeffs = (int *)malloc((p1->degree + 1) * sizeof(int));
		assert(coeffs);
		int i = 0;
		for (; i <= p2->degree; i++)
			coeffs[i] = (p1->coeffs[i] >= p2->coeffs[i]) ? (p1->coeffs[i] - p2->coeffs[i]) : (p1->coeffs[i] - p2->coeffs[i] + p);
		for (; i <= p1->degree; i++)
			coeffs[i] = p1->coeffs[i];
		set_poly(dif, p1->degree, coeffs);
	}
	else if (p1->degree < p2->degree)
	{
		int *coeffs = (int *)malloc((p2->degree + 1) * sizeof(int));
		assert(coeffs);
		int i = 0;
		for (; i <= p1->degree; i++)
			coeffs[i] = (p1->coeffs[i] >= p2->coeffs[i]) ? (p1->coeffs[i] - p2->coeffs[i]) : (p1->coeffs[i] - p2->coeffs[i] + p);
		for (; i <= p2->degree; i++)
			coeffs[i] = (p2->coeffs[i] == 0) ? (0) : (p - p2->coeffs[i]);
		set_poly(dif, p2->degree, coeffs);
	}
	else if (p1->degree == -1)
	{
		int *coeff = (int *)calloc(1, sizeof(int));
		assert(coeff);
		set_poly(dif, -1, coeff);
	}
	else
	{
		int i;
		for (i = p1->degree; i >= 0 && (p1->coeffs[i] - p2->coeffs[i]) % p == 0; i--)
			;
		if (i == -1)
		{
			int *coeff = (int *)calloc(1, sizeof(int));
			assert(coeff);
			set_poly(dif, -1, coeff);
		}
		else
		{
			int degree = i;
			int *coeffs = (int *)malloc((degree + 1) * sizeof(int));
			assert(coeffs);
			for (; i >= 0; i--)
				coeffs[i] = (p1->coeffs[i] >= p2->coeffs[i]) ? (p1->coeffs[i] - p2->coeffs[i]) : (p1->coeffs[i] - p2->coeffs[i] + p);
			set_poly(dif, degree, coeffs);
		}
	}
	return dif;
}

poly *mul_poly(poly *p1, poly *p2)
{
	poly *prod = new_poly();
	if (p1->degree == -1 || p2->degree == -1)
	{
		int *coeff = (int *)calloc(1, sizeof(int));
		assert(coeff);
		set_poly(prod, -1, coeff);
	}
	else
	{
		int degree = p1->degree + p2->degree;
		int *coeffs = (int *)calloc(degree + 1, sizeof(int));
		assert(coeffs);
		for (int i = 0; i <= p1->degree; i++)
			for (int j = 0; j <= p2->degree; j++)
				coeffs[i + j] = (coeffs[i + j] + (p1->coeffs[i] * p2->coeffs[j])) % p;
		set_poly(prod, degree, coeffs);
	}
	return prod;
}

int inverse_zp(int n)
{
	int v0 = 0;
	int r0 = p;
	int v1 = 1;
	int r1 = n;
	while (r1 != 1)
	{
		int q = r0 / r1;
		int v2 = v0 - q * v1;
		int r2 = r0 - q * r1;
		v0 = v1;
		r0 = r1;
		v1 = v2;
		r1 = r2;
	}
	v1 %= p;
	if (v1 < 0)
		v1 += p;
	return v1;
}

void euclid_division(poly *p1, poly *p2, poly *q, poly *r)
{
	int *rest = (int *)malloc((p1->degree + 1) * sizeof(int));
	assert(rest);
	for (int i = 0; i <= p1->degree; i++)
		rest[i] = p1->coeffs[i];
	q->degree = p1->degree - p2->degree;
	q->coeffs = (int *)calloc(q->degree + 1, sizeof(int));
	assert(q->coeffs);
	int rest_degree = p1->degree;
	int inv_LC_p2 = inverse_zp(p2->coeffs[p2->degree]);
	while (rest_degree >= p2->degree)
	{
		int coeff = (rest[rest_degree] * inv_LC_p2) % p;
		q->coeffs[rest_degree - p2->degree] = coeff;
		for (int i = 0; i <= p2->degree; i++)
		{
			rest[i + rest_degree - p2->degree] -= coeff * p2->coeffs[i];
			rest[i + rest_degree - p2->degree] %= p;
			if (rest[i + rest_degree - p2->degree] < 0)
				rest[i + rest_degree - p2->degree] += p;
		}
		while (rest_degree >= 0 && rest[rest_degree] == 0)
			rest_degree--;
	}
	r->degree = rest_degree;
	if (r->degree == -1)
		r->coeffs = (int *)calloc(1, sizeof(int));
	else
		r->coeffs = (int *)malloc((r->degree + 1) * sizeof(int));
	assert(r->coeffs);
	for (int i = 0; i <= r->degree; i++)
		r->coeffs[i] = rest[i];
	free(rest);
}

poly *xgcd(poly *p1, poly *p2, poly *u, poly *v)
{
	// ui * p1 + vi * p2 = ri
	poly *u0 = new_poly_1();
	poly *v0 = new_poly_0();
	poly *r0 = new_poly_from_copy(p1);
	poly *u1 = new_poly_0();
	poly *v1 = new_poly_1();
	poly *r1 = new_poly_from_copy(p2);
	while (r1->degree != -1)
	{
		poly *q = new_poly();
		poly *r2 = new_poly();
		euclid_division(r0, r1, q, r2);
		poly *qu1 = mul_poly(q, u1);
		poly *u2 = substract_poly(u0, qu1);
		poly *qv1 = mul_poly(q, v1);
		poly *v2 = substract_poly(v0, qv1);
		free_poly_full(q);
		free_poly_full(qu1);
		free_poly_full(qv1);
		free_poly_full(u0);
		free_poly_full(v0);
		free_poly_full(r0);
		u0 = u1;
		v0 = v1;
		r0 = r1;
		u1 = u2;
		v1 = v2;
		r1 = r2;
	}
	u = u0;
	v = v0;
	free_poly_full(u1);
	free_poly_full(v1);
	free_poly_full(r1);
	return r0;
}



