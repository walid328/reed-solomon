#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "finite_field.h"

int p = 5;


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

void str_add_int(char* string, int *index, int n)
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

char *string_polynomial_iso_length(poly *q)
{
	if (q->degree == -1)
	{
		char *poly_string = (char*)malloc(2 * sizeof(char));
		assert(poly_string != NULL);
		poly_string[0] = '0';
		poly_string[1] = '\0';
		return poly_string;
	}
	int size_p = size_int_str(p);
	int size_d = size_int_str(q->degree);
	char *poly_string = (char*)malloc(((size_p + size_d + 6) * (q->degree + 1) + 1) * sizeof(char));
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


char *string_polynomial_minimal(poly *q)
{
	if (q->degree == -1)
	{
		char *poly_string = (char*)malloc(2 * sizeof(char));
		assert(poly_string != NULL);
		poly_string[0] = '0';
		poly_string[1] = '\0';
		return poly_string;
	}
	int size_p = size_int_str(p);
	int size_d = size_int_str(q->degree);
	char *poly_string = (char*)malloc(((size_p + size_d + 6) * (q->degree + 1) + 1) * sizeof(char));
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

void print_polynomial_iso_length(poly *q)
{
	char *poly_string = string_polynomial_iso_length(q);
    printf("%s\n", poly_string);
	free(poly_string);
}

void print_polynomial_minimal(poly *q)
{
	char *poly_string = string_polynomial_minimal(q);
    printf("%s\n", poly_string);
	free(poly_string);
}

void print_polynomial(poly *q)
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

poly *add_polynomials(poly *p1, poly *p2)
{
    poly *sum = (poly *)malloc(sizeof(poly));
    assert(sum != NULL);
    if (p1->degree > p2->degree)
    {
        sum->degree = p1->degree;
        sum->coefficients = (int *)malloc(sum->degree * sizeof(int));
        assert(sum->coefficients != NULL);
        int i = 0;
        for (; i <= p2->degree; i++)
            sum->coefficients[i] = (p1->coefficients[i] + p2->coefficients[i]) % p;
        for (; i <= p1->degree; i++)
            sum->coefficients[i] = p1->coefficients[i];
    }
    else if (p1->degree < p2->degree)
    {
        sum->degree = p2->degree;
        sum->coefficients = (int *)malloc(sum->degree * sizeof(int));
        assert(sum->coefficients != NULL);
        int i = 0;
        for (; i <= p1->degree; i++)
            sum->coefficients[i] = (p1->coefficients[i] + p2->coefficients[i]) % p;
        for (; i <= p2->degree; i++)
            sum->coefficients[i] = p2->coefficients[i];
    }
    else if (p1->degree == -1)
    {
        sum->degree = -1;
        sum->coefficients = NULL;
    }
    else
    {
        int i;
        for (i = p1->degree; i >= 0 && (p1->coefficients[i] + p2->coefficients[i]) % p == 0; i--)
            ;
        sum->degree = i;
        if (i == -1)
            sum->coefficients = NULL;
        else
        {
            sum->coefficients = (int *)malloc(sum->degree * sizeof(int));
            assert(sum->coefficients != NULL);
            for (; i >= 0; i--)
                sum->coefficients[i] = (p1->coefficients[i] + p2->coefficients[i]) % p;
        }
    }
    return sum;
}