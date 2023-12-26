#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "rs_code.h"
#include "polynomial.h"
#include "finite_field.h"

void encoding(int *c, int *a, int n, int *m, int k)
{
    poly *f = new_poly();
    set_poly(f, k - 1, m);
    multi_eval_poly(f, n, a, c);
    free_poly(f);
}

void decoding(int *m, int *a, int *b, int n, int k)
{
    poly *g_0 = new_poly_1();
    for (int i = 0; i < n; i++)
    {
        poly *x_m_ai = new_poly_from_coeffs(1, opp_zp(a[i]), 1);
        mul_poly(g_0, g_0, x_m_ai);
        free_full_poly(x_m_ai);
    }
    poly *g_1 = interpolation(a, b, n);
    poly *g, *u, *v;
    partial_gcd_poly(&g, &u, &v, g_0, g_1, (n + k) / 2);
    poly *f_1 = new_poly();
    poly *r = new_poly();
    euclid_div_poly(f_1, r, g, v);
    if (degree(r) == -1 && degree(f_1) < k)
    {
        for (int i = 0; i < k; i++)
            m[i] = coeffs(f_1)[i];
        free_full_poly(g_0);
        free_full_poly(g_1);
        free_full_poly(g);
        free_full_poly(u);
        free_full_poly(v);
        free_full_poly(r);
        free_full_poly(f_1);
    }
    else
    {
        free_full_poly(g_0);
        free_full_poly(g_1);
        free_full_poly(g);
        free_full_poly(u);
        free_full_poly(v);
        free_full_poly(r);
        free_full_poly(f_1);
        fprintf(stderr, "Decoding failure: too many errors!\n");
        exit(EXIT_FAILURE);
    }
}