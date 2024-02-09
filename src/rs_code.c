#include <stdlib.h>
#include <stdio.h>

#include "rs_code.h"
#include "finite_field.h"
#include "field_settings.h"
#include "zp_array.h"
#include "polynomial.h"

/******************************************************/

void rs_g_0(poly g_0, zp_array points, int block_length)
{
    poly_reset_coeffs(g_0, 0, 1);
    for (int i = 0; i < block_length; i++)
    {
        poly x_m_ai = poly_new_set(1, zp_opp(points[i]), 1);
        poly_mul(g_0, g_0, x_m_ai);
        poly_free(x_m_ai);
    }
}

void rs_g_0_fourier(poly g_0, int block_length)
{
    poly_reset_coeffs(g_0, 0, 1);
    for (int i = 0; i < block_length; i++)
    {
        poly x_m_ai = poly_new_set(1, zp_opp(omegas[i * n / block_length]), 1);
        poly_mul(g_0, g_0, x_m_ai);
        poly_free(x_m_ai);
    }
}

/******************************************************/

zp_array rs_encode(int block_length, int message_length, zp_array points, zp_array message)
{
    poly f = poly_new();
    poly_set(f, message_length - 1, message);
    zp_array codeword = poly_eval_zp_array(f, points, block_length);
    poly_set(f, -1, NULL);
    poly_free(f);
    return codeword;
}

zp_array rs_encode_2(int block_length, int message_length, zp_array message)
{
    poly f = poly_new();
    poly_set(f, message_length - 1, message);
    zp_array codeword = poly_dft(f, block_length);
    poly_set(f, -1, NULL);
    poly_free(f);
    return codeword;
}

zp_array rs_fast_encode(int block_length, int message_length, zp_array message)
{
    poly f = poly_new();
    poly_set(f, message_length - 1, message);
    zp_array codeword = poly_fft(f, block_length);
    poly_set(f, -1, NULL);
    poly_free(f);
    return codeword;
}

/******************************************************/

zp_array rs_decode(poly g_0, int block_length, int message_length, zp_array points, zp_array received)
{
    poly g_1 = poly_new();
    interpolation(g_1, points, received, block_length);
    poly g = poly_new();
    poly u = poly_new();
    poly v = poly_new();
    poly_xgcd_partial(g, u, v, g_0, g_1, (block_length + message_length + 1) / 2);
    poly f_1 = poly_new();
    poly r = poly_new();
    poly_euc_div(f_1, r, g, v);
    poly_free_multi(4, g_1, g, u, v);
    if (poly_is_zero(r) && poly_deg(f_1) < message_length)
    {
        zp_array message = zp_array_new_zeros(message_length);
        for (int i = 0; i <= poly_deg(f_1); i++)
            message[i] = poly_coeffs(f_1)[i];
        poly_free_multi(2, r, f_1);
        return message;
    }
    else
    {
        poly_free_multi(2, r, f_1);
        fprintf(stderr, "rs_decode failure: too many errors!\n");
        exit(EXIT_FAILURE);
    }
}

zp_array rs_decode_2(poly g_0, int block_length, int message_length, zp_array received)
{
    poly g_1 = poly_new();
    poly_inv_dft(g_1, received, block_length);
    poly g = poly_new();
    poly u = poly_new();
    poly v = poly_new();
    poly_xgcd_partial(g, u, v, g_0, g_1, (block_length + message_length + 1) / 2);
    poly f_1 = poly_new();
    poly r = poly_new();
    poly_euc_div(f_1, r, g, v);
    poly_free_multi(4, g_1, g, u, v);
    if (poly_is_zero(r) && poly_deg(f_1) < message_length)
    {
        zp_array message = zp_array_new_zeros(message_length);
        for (int i = 0; i <= poly_deg(f_1); i++)
            message[i] = poly_coeffs(f_1)[i];
        poly_free_multi(2, r, f_1);
        return message;
    }
    else
    {
        poly_free_multi(2, r, f_1);
        fprintf(stderr, "rs_decode_2 failure: too many errors!\n");
        exit(EXIT_FAILURE);
    }
}

zp_array rs_fast_decode(poly g_0, int block_length, int message_length, zp_array received)
{
    poly g_1 = poly_new();
    poly_inv_fft(g_1, received, block_length);
    poly g = poly_new();
    poly u = poly_new();
    poly v = poly_new();
    poly_fast_xgcd_partial(g, u, v, g_0, g_1, (block_length + message_length + 1) / 2);
    poly f_1 = poly_new();
    poly r = poly_new();
    poly_fast_euc_div(f_1, r, g, v);
    poly_free_multi(4, g_1, g, u, v);
    if (poly_is_zero(r) && poly_deg(f_1) < message_length)
    {
        zp_array message = zp_array_new_zeros(message_length);
        for (int i = 0; i <= poly_deg(f_1); i++)
            message[i] = poly_coeffs(f_1)[i];
        poly_free_multi(2, r, f_1);
        return message;
    }
    else
    {
        poly_free_multi(2, r, f_1);
        fprintf(stderr, "rs_fast_decode failure: too many errors!\n");
        exit(EXIT_FAILURE);
    }
}