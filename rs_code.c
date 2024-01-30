#include <stdlib.h>
#include <stdio.h>

#include "rs_code.h"
#include "finite_field.h"
#include "field_settings.h"
#include "array.h"
#include "polynomial.h"

array rs_encode(int n, int k, array points, array message)
{
    poly f = poly_new();
    poly_set(f, k - 1, message);
    array codeword = poly_eval_array(f, points, n);
    poly_set(f, -1, NULL);
    poly_free(f);
    return codeword;
}

array rs_decode(int n, int k, array points, array received)
{
    poly g_0 = poly_new_set(0, 1);
    for (int i = 0; i < n; i++)
    {
        poly x_m_ai = poly_new_set(1, zp_opp(points[i]), 1);
        poly_mul(g_0, g_0, x_m_ai);
        poly_free(x_m_ai);
    }
    poly g_1 = poly_new();
    interpolation(g_1, points, received, n);
    poly g = poly_new();
    poly u = poly_new();
    poly v = poly_new();
    poly_xgcd_partial(&g, &u, &v, g_0, g_1, (n + k) / 2);
    poly f_1 = poly_new();
    poly r = poly_new();
    poly_euc_div(f_1, r, g, v);
    poly_free(g_0);
    poly_free(g_1);
    poly_free(g);
    poly_free(u);
    poly_free(v);
    if (poly_deg(r) == -1 && poly_deg(f_1) < k)
    {
        array message = array_new(k);
        for (int i = 0; i < k; i++)
            message[i] = poly_coeffs(f_1)[i];
        poly_free(r);
        poly_free(f_1);
        return message;
    }
    else
    {
        poly_free(r);
        poly_free(f_1);
        fprintf(stderr, "rs_decode failure: too many errors!\n");
        exit(EXIT_FAILURE);
    }
}