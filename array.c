#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

#include "array.h"
#include "field_settings.h"
#include "finite_field.h"

/******************************************************/

void array_print(array tab, int tab_size)
{
    for (int i = 0; i < tab_size; i++)
        printf("%d ", tab[i]);
    printf("\n");
}

/******************************************************/

array array_new(int tab_size)
{
    array tab = (array)malloc(tab_size * sizeof(int));
    assert(tab);
    return tab;
}

array array_new_zeros(int tab_size)
{
    array tab = (array)calloc(tab_size, sizeof(int));
    assert(tab);
    return tab;
}

array array_new_set(int tab_size, ...)
{
    array tab = array_new(tab_size);
    va_list valist;
    va_start(valist, tab_size);
    for (int i = 0; i < tab_size; i++)
        tab[i] = va_arg(valist, int);
    va_end(valist);
    return tab;
}

array array_new_rand(int tab_size)
{
    array tab = array_new(tab_size);
    for (int i = 0; i < tab_size; i++)
        tab[i] = zp_rand();
    return tab;
}

/******************************************************/

void array_free(array tab)
{
    free(tab);
}

/******************************************************/

void array_set(array tab, int tab_size, ...)
{
    va_list valist;
    va_start(valist, tab_size);
    for (int i = 0; i < tab_size; i++)
        tab[i] = va_arg(valist, int);
    va_end(valist);
}

bool array_equal(array tab1, array tab2, int tab_size)
{
    for (int i = 0; i < tab_size; i++)
        if (tab1[i] != tab2[i])
            return false;
    return true;
}

void array_add_errors(array tab, int tab_size, int number_errors)
{
    for (int i = 0; i < number_errors; i++)
    {
        int random_position = rand() % tab_size;
        tab[random_position] = zp_rand();
    }
}

/******************************************************/

void array_split(array *even, array *odd, array tab, int tab_size)
{
    *even = array_new(tab_size / 2);
    *odd = array_new(tab_size / 2);
    for (int i = 0; i < tab_size / 2; i++)
    {
        (*even)[i] = tab[2 * i];
        (*odd)[i] = tab[2 * i + 1];
    }
}

array array_merge(array even, array odd, int subtab_size)
{
    array tab = array_new(2 * subtab_size);
    for (int i = 0; i < subtab_size; i++)
    {
        tab[2 * i] = even[i];
        tab[2 * i + 1] = odd[i];
    }
    return tab;
}

/******************************************************/

array array_fft(array coeffs, int d)
{
    if (d == 1)
    {
        array eval = array_new(1);
        eval[0] = coeffs[0];
        return eval;
    }
    int half_d = d / 2;
    array r0 = array_new(half_d);
    array r1 = array_new(half_d);
    for (int i = 0; i < half_d; i++)
    {
        r0[i] = zp_add(coeffs[i], coeffs[i + half_d]);
        r1[i] = zp_mul(zp_sub(coeffs[i], coeffs[i + half_d]), omegas[i * (n / d)]);
    }
    array eval_r0 = array_fft(r0, half_d);
    array eval_r1 = array_fft(r1, half_d);
    array_free(r0);
    array_free(r1);
    array eval = array_merge(eval_r0, eval_r1, half_d);
    array_free(eval_r0);
    array_free(eval_r1);
    return eval;
}

array array_inv_fft(array eval, int d)
{
    if (d == 1)
    {
        array coeffs = array_new(1);
        coeffs[0] = eval[0];
        return coeffs;
    }
    int half_d = d / 2;
    array eval_r0, eval_r1;
    array_split(&eval_r0, &eval_r1, eval, d);
    array r0 = array_inv_fft(eval_r0, half_d);
    array r1 = array_inv_fft(eval_r1, half_d);
    array_free(eval_r0);
    array_free(eval_r1);
    for (int i = 1; i <= half_d; i++)
        r1[half_d - i] = zp_opp(zp_mul(r1[half_d - i], omegas[i * (n / d)]));
    array coeffs = array_new(d);
    for (int i = 0; i < half_d; i++)
    {
        coeffs[i] = zp_mul(zp_add(r0[i], r1[i]), (p + 1) / 2);
        coeffs[i + half_d] = zp_mul(zp_sub(r0[i], r1[i]), (p + 1) / 2);
    }
    array_free(r0);
    array_free(r1);
    return coeffs;
}

/******************************************************/

array formal_serie_mul(array P, array Q, int d)
{
	if (d == 2 * n)
	{
		array P0 = array_new_zeros(n);
		array P1 = array_new_zeros(n);
		array Q0 = array_new_zeros(n);
		array Q1 = array_new_zeros(n);
		for (int i = 0; i < n / 2; i++)
			P0[i] = P[i];
		for (int i = n / 2; i < n; i++)
			P1[i - n / 2] = P[i];
		for (int i = 0; i < n / 2; i++)
			Q0[i] = Q[i];
		for (int i = n / 2; i < n; i++)
			Q1[i - n / 2] = Q[i];
		array P0Q0 = formal_serie_mul(P0, Q0, n);
		array P0Q1 = formal_serie_mul(P0, Q1, n);
		array P1Q0 = formal_serie_mul(P1, Q0, n);
		array res = array_new_zeros(d);
		for (int i = 0; i < n / 2; i++)
			res[i] = P0Q0[i];
		for (int i = n / 2; i < n; i++)
			res[i] = zp_add(P0Q0[i], zp_add(P0Q1[i - n / 2], P1Q0[i - n / 2]));
		array_free(P0);
		array_free(P1);
		array_free(Q0);
		array_free(Q1);
		array_free(P0Q0);
		array_free(P0Q1);
		array_free(P1Q0);
		return res;
	}
    array eval_P = array_fft(P, d);
    array eval_Q = array_fft(Q, d);
    array eval_res = array_new(d);
    for (int i = 0; i < d; i++)
        eval_res[i] = zp_mul(eval_P[i], eval_Q[i]);
    array res = array_inv_fft(eval_res, d);
    array_free(eval_P);
    array_free(eval_Q);
    array_free(eval_res);
    return res;
}

array formal_serie_inv(array P, int d)
{
    if (d == 1)
    {
        array Q = array_new(1);
        Q[0] = zp_inv(P[0]);
        return Q;
    }
    array Q0 = formal_serie_inv(P, d / 2);
    array Q0_d = array_new_zeros(d);
    for (int i = 0; i < d / 2; i++)
        Q0_d[i] = Q0[i];
    array Qa = formal_serie_mul(Q0_d, Q0_d, d);
    array P_2d = array_new_zeros(2 * d);
    for (int i = 0; i < d; i++)
        P_2d[i] = P[i];
	array Qa_2d = array_new_zeros(2 * d);
	for (int i = 0; i < d; i++)
		Qa_2d[i] = Qa[i];
    array Qa2 = formal_serie_mul(P_2d, Qa_2d, 2 * d);
    array Q = array_new_zeros(d);
    for (int i = 0; i < d / 2; i++)
        Q[i] = Q0[i];
    for (int i = d / 2; i < d; i++)
        Q[i] = zp_opp(Qa2[i]);
    array_free(Q0);
    array_free(Q0_d);
    array_free(Qa);
    array_free(P_2d);
    array_free(Qa_2d);
    array_free(Qa2);
    return Q;
}

void formal_serie_fast_euc_div(array *quo, array *rem, array P, int deg_p, array D, int deg_d)
{
    // dm is the smallest power of 2 greater or equal to deg_d
    int dm = 1;
    while (dm < deg_d)
        dm <<= 1;
    int dnmm = 1;
    // dnmm is the smallest power of 2 greater or equal to op_p - deg_d + 1
    while (dnmm < deg_p - deg_d + 1)
        dnmm <<= 1;
    array P_prime = array_new_zeros(2 * dnmm);
    array D_prime = array_new_zeros(dnmm);
    for (int i = 0; i < deg_p - deg_d + 1; i++)
        P_prime[i] = P[deg_p - i];
    for (int i = 0; (i < (deg_p - deg_d + 1)) && (i < (deg_d + 1)); i++)
        D_prime[i] = D[deg_d - i];
    array D_prime_prime = formal_serie_inv(D_prime, dnmm);
	array D_prime_prime_2dnmm = array_new_zeros(2 * dnmm);
	for (int i = 0; i < dnmm; i++)
		D_prime_prime_2dnmm[i] = D_prime_prime[i];
    array Q_prime = formal_serie_mul(P_prime, D_prime_prime_2dnmm, 2 * dnmm);
    *quo = array_new(deg_p - deg_d + 1);
    for (int i = 0; i < deg_p - deg_d + 1; i++)
        (*quo)[i] = Q_prime[deg_p - deg_d - i];
    array Q_2dm = array_new_zeros(2 * dm);
    for (int i = 0; i < deg_p - deg_d + 1 && i < deg_d; i++)
        Q_2dm[i] = (*quo)[i];
    array D_2dm = array_new_zeros(2 * dm);
    for (int i = 0; i < deg_d; i++)
        D_2dm[i] = D[i];
    array QD = formal_serie_mul(Q_2dm, D_2dm, 2 * dm);
    *rem = array_new(deg_d);
    for (int i = 0; i < deg_d; i++)
        (*rem)[i] = zp_sub(P[i], QD[i]);
    array_free(P_prime);
    array_free(D_prime);
    array_free(D_prime_prime);
    array_free(D_prime_prime_2dnmm);
    array_free(Q_prime);
    array_free(Q_2dm);
    array_free(D_2dm);
    array_free(QD);
}