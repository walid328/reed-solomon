#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

#include "zp_array.h"
#include "field_settings.h"

/******************************************************/

void zp_array_print(zp_array tab, int tab_size)
{
    for (int i = 0; i < tab_size; i++)
        printf("%d ", tab[i]);
    printf("\n");
}

/******************************************************/

zp_array zp_array_new(int tab_size)
{
    zp_array tab = (zp_array)malloc(tab_size * sizeof(zp_t));
    assert(tab);
    return tab;
}

zp_array zp_array_new_zeros(int tab_size)
{
    zp_array tab = (zp_array)calloc(tab_size, sizeof(zp_t));
    assert(tab);
    return tab;
}

zp_array zp_array_new_set(int tab_size, ...)
{
    zp_array tab = zp_array_new(tab_size);
    va_list valist;
    va_start(valist, tab_size);
    for (int i = 0; i < tab_size; i++)
        tab[i] = zp_mod(va_arg(valist, int));
    va_end(valist);
    return tab;
}

zp_array zp_array_new_rand(int tab_size)
{
    zp_array tab = zp_array_new(tab_size);
    for (int i = 0; i < tab_size; i++)
        tab[i] = zp_rand();
    return tab;
}

/******************************************************/

void zp_array_free(zp_array tab)
{
    free(tab);
}

/******************************************************/

void zp_array_set(zp_array tab, int tab_size, ...)
{
    va_list valist;
    va_start(valist, tab_size);
    for (int i = 0; i < tab_size; i++)
        tab[i] = zp_mod(va_arg(valist, int));
    va_end(valist);
}

bool zp_array_equal(zp_array tab1, zp_array tab2, int tab_size)
{
    for (int i = 0; i < tab_size; i++)
        if (tab1[i] != tab2[i])
            return false;
    return true;
}

void zp_array_add_errors(zp_array tab, int tab_size, int number_errors)
{
    for (int i = 0; i < number_errors; i++)
    {
        int random_position = rand() % tab_size;
        tab[random_position] = zp_rand();
    }
}

/******************************************************/

/* Split and merge */

// Split an zp_array with even an even number of int
// by even and odd indexes.
void zp_array_split(zp_array *even, zp_array *odd, zp_array tab, int tab_size)
{
    *even = zp_array_new(tab_size / 2);
    *odd = zp_array_new(tab_size / 2);
    for (int i = 0; i < tab_size / 2; i++)
    {
        (*even)[i] = tab[2 * i];
        (*odd)[i] = tab[2 * i + 1];
    }
}

// Merge two zp_arrays into one. Even contains int of
// even indexes and odd containts int of odd indexes.
zp_array zp_array_merge(zp_array even, zp_array odd, int subtab_size)
{
    zp_array tab = zp_array_new(2 * subtab_size);
    for (int i = 0; i < subtab_size; i++)
    {
        tab[2 * i] = even[i];
        tab[2 * i + 1] = odd[i];
    }
    return tab;
}

/******************************************************/

zp_array zp_array_fft(zp_array coeffs, int d)
{
    if (d == 1)
    {
        zp_array eval = zp_array_new(1);
        eval[0] = coeffs[0];
        return eval;
    }
    int half_d = d / 2;
    zp_array r0 = zp_array_new(half_d);
    zp_array r1 = zp_array_new(half_d);
    for (int i = 0; i < half_d; i++)
    {
        r0[i] = zp_add(coeffs[i], coeffs[i + half_d]);
        r1[i] = zp_mul(zp_sub(coeffs[i], coeffs[i + half_d]), omegas[i * (n / d)]);
    }
    zp_array eval_r0 = zp_array_fft(r0, half_d);
    zp_array eval_r1 = zp_array_fft(r1, half_d);
    zp_array_free(r0);
    zp_array_free(r1);
    zp_array eval = zp_array_merge(eval_r0, eval_r1, half_d);
    zp_array_free(eval_r0);
    zp_array_free(eval_r1);
    return eval;
}

zp_array zp_array_inv_fft(zp_array eval, int d)
{
    if (d == 1)
    {
        zp_array coeffs = zp_array_new(1);
        coeffs[0] = eval[0];
        return coeffs;
    }
    int half_d = d / 2;
    zp_array eval_r0, eval_r1;
    zp_array_split(&eval_r0, &eval_r1, eval, d);
    zp_array r0 = zp_array_inv_fft(eval_r0, half_d);
    zp_array r1 = zp_array_inv_fft(eval_r1, half_d);
    zp_array_free(eval_r0);
    zp_array_free(eval_r1);
    for (int i = 1; i <= half_d; i++)
        r1[half_d - i] = zp_opp(zp_mul(r1[half_d - i], omegas[i * (n / d)]));
    zp_array coeffs = zp_array_new(d);
    for (int i = 0; i < half_d; i++)
    {
        coeffs[i] = zp_mul(zp_add(r0[i], r1[i]), (p + 1) / 2);
        coeffs[i + half_d] = zp_mul(zp_sub(r0[i], r1[i]), (p + 1) / 2);
    }
    zp_array_free(r0);
    zp_array_free(r1);
    return coeffs;
}

/******************************************************/

zp_array formal_serie_mul(zp_array P, zp_array Q, int d)
{
    if (d == 2 * n)
    {
        zp_array P0 = zp_array_new_zeros(n);
        zp_array P1 = zp_array_new_zeros(n);
        zp_array Q0 = zp_array_new_zeros(n);
        zp_array Q1 = zp_array_new_zeros(n);
        for (int i = 0; i < n / 2; i++)
            P0[i] = P[i];
        for (int i = n / 2; i < n; i++)
            P1[i - n / 2] = P[i];
        for (int i = 0; i < n / 2; i++)
            Q0[i] = Q[i];
        for (int i = n / 2; i < n; i++)
            Q1[i - n / 2] = Q[i];
        zp_array P0Q0 = formal_serie_mul(P0, Q0, n);
        zp_array P0Q1 = formal_serie_mul(P0, Q1, n);
        zp_array P1Q0 = formal_serie_mul(P1, Q0, n);
        zp_array res = zp_array_new_zeros(d);
        for (int i = 0; i < n / 2; i++)
            res[i] = P0Q0[i];
        for (int i = n / 2; i < n; i++)
            res[i] = zp_add(P0Q0[i], zp_add(P0Q1[i - n / 2], P1Q0[i - n / 2]));
        zp_array_free(P0);
        zp_array_free(P1);
        zp_array_free(Q0);
        zp_array_free(Q1);
        zp_array_free(P0Q0);
        zp_array_free(P0Q1);
        zp_array_free(P1Q0);
        return res;
    }
    zp_array eval_P = zp_array_fft(P, d);
    zp_array eval_Q = zp_array_fft(Q, d);
    zp_array eval_res = zp_array_new(d);
    for (int i = 0; i < d; i++)
        eval_res[i] = zp_mul(eval_P[i], eval_Q[i]);
    zp_array res = zp_array_inv_fft(eval_res, d);
    zp_array_free(eval_P);
    zp_array_free(eval_Q);
    zp_array_free(eval_res);
    return res;
}

zp_array formal_serie_inv(zp_array P, int d)
{
    if (d == 1)
    {
        zp_array Q = zp_array_new(1);
        Q[0] = zp_inv(P[0]);
        return Q;
    }
    zp_array Q0 = formal_serie_inv(P, d / 2);
    zp_array Q0_d = zp_array_new_zeros(d);
    for (int i = 0; i < d / 2; i++)
        Q0_d[i] = Q0[i];
    zp_array Qa = formal_serie_mul(Q0_d, Q0_d, d);
    zp_array P_2d = zp_array_new_zeros(2 * d);
    for (int i = 0; i < d; i++)
        P_2d[i] = P[i];
    zp_array Qa_2d = zp_array_new_zeros(2 * d);
    for (int i = 0; i < d; i++)
        Qa_2d[i] = Qa[i];
    zp_array Qa2 = formal_serie_mul(P_2d, Qa_2d, 2 * d);
    zp_array Q = zp_array_new_zeros(d);
    for (int i = 0; i < d / 2; i++)
        Q[i] = Q0[i];
    for (int i = d / 2; i < d; i++)
        Q[i] = zp_opp(Qa2[i]);
    zp_array_free(Q0);
    zp_array_free(Q0_d);
    zp_array_free(Qa);
    zp_array_free(P_2d);
    zp_array_free(Qa_2d);
    zp_array_free(Qa2);
    return Q;
}

void formal_serie_fast_euc_div(zp_array *quo, zp_array *rem, zp_array P, int deg_p, zp_array D, int deg_d)
{
    // dm is the smallest power of 2 greater or equal to deg_d
    int dm = 1;
    while (dm < deg_d)
        dm <<= 1;
    int dnmm = 1;
    // dnmm is the smallest power of 2 greater or equal to op_p - deg_d + 1
    while (dnmm < deg_p - deg_d + 1)
        dnmm <<= 1;
    zp_array P_prime = zp_array_new_zeros(2 * dnmm);
    zp_array D_prime = zp_array_new_zeros(dnmm);
    for (int i = 0; i < deg_p - deg_d + 1; i++)
        P_prime[i] = P[deg_p - i];
    for (int i = 0; (i < (deg_p - deg_d + 1)) && (i < (deg_d + 1)); i++)
        D_prime[i] = D[deg_d - i];
    zp_array D_prime_prime = formal_serie_inv(D_prime, dnmm);
    zp_array D_prime_prime_2dnmm = zp_array_new_zeros(2 * dnmm);
    for (int i = 0; i < dnmm; i++)
        D_prime_prime_2dnmm[i] = D_prime_prime[i];
    zp_array Q_prime = formal_serie_mul(P_prime, D_prime_prime_2dnmm, 2 * dnmm);
    *quo = zp_array_new(deg_p - deg_d + 1);
    for (int i = 0; i < deg_p - deg_d + 1; i++)
        (*quo)[i] = Q_prime[deg_p - deg_d - i];
    zp_array Q_2dm = zp_array_new_zeros(2 * dm);
    for (int i = 0; i < deg_p - deg_d + 1 && i < deg_d; i++)
        Q_2dm[i] = (*quo)[i];
    zp_array D_2dm = zp_array_new_zeros(2 * dm);
    for (int i = 0; i < deg_d; i++)
        D_2dm[i] = D[i];
    zp_array QD = formal_serie_mul(Q_2dm, D_2dm, 2 * dm);
    *rem = zp_array_new(deg_d);
    for (int i = 0; i < deg_d; i++)
        (*rem)[i] = zp_sub(P[i], QD[i]);
    zp_array_free(P_prime);
    zp_array_free(D_prime);
    zp_array_free(D_prime_prime);
    zp_array_free(D_prime_prime_2dnmm);
    zp_array_free(Q_prime);
    zp_array_free(Q_2dm);
    zp_array_free(D_2dm);
    zp_array_free(QD);
}