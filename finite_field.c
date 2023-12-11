#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "finite_field.h"

int p = 5;

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