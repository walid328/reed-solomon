#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "rs_code.h"
#include "polynomial.h"
#include "finite_field.h"

void rs_encode(int *c, int *a, int n, int *m, int k)
{
    poly *f = poly_new();
    poly_set(f, k - 1, m);
    poly_eval_multi(f, n, a, c);
    poly_free(f);
}

void rs_decode(int *m, int *a, int *b, int n, int k)
{
    poly *g_0 = poly_new_1();
    for (int i = 0; i < n; i++)
    {
        poly *x_m_ai = poly_new_from_poly_coeffs(1, zp_opp(a[i]), 1);
        poly_mul(g_0, g_0, x_m_ai);
        poly_free_full(x_m_ai);
    }
    poly *g_1 = interpolation(a, b, n);
    poly *g, *u, *v;
    poly_xgcd_partial(&g, &u, &v, g_0, g_1, (n + k) / 2);
    poly *f_1 = poly_new();
    poly *r = poly_new();
    poly_euc_div(f_1, r, g, v);
    if (poly_degree(r) == -1 && poly_degree(f_1) < k)
    {
        for (int i = 0; i < k; i++)
            m[i] = poly_coeffs(f_1)[i];
        poly_free_full(g_0);
        poly_free_full(g_1);
        poly_free_full(g);
        poly_free_full(u);
        poly_free_full(v);
        poly_free_full(r);
        poly_free_full(f_1);
    }
    else
    {
        poly_free_full(g_0);
        poly_free_full(g_1);
        poly_free_full(g);
        poly_free_full(u);
        poly_free_full(v);
        poly_free_full(r);
        poly_free_full(f_1);
        fprintf(stderr, "rs_decode failure: too many errors!\n");
        exit(EXIT_FAILURE);
    }
}