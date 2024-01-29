#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>

#include "polynomial.h"
#include "field_settings.h"
#include "finite_field.h"
#include "array.h"

/******************************************************/

struct polynomial
{
    int deg;
    array coeffs;
};

/******************************************************/

int poly_deg(const poly f)
{
    return f->deg;
}

array poly_coeffs(const poly f)
{
    return f->coeffs;
}

void poly_set(poly f, int deg, array coeffs)
{
    f->deg = deg;
    f->coeffs = coeffs;
}

/******************************************************/

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

void str_add_int(char *string, array index, int n)
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

char *str_poly_iso_length(poly f)
{
    if (f->deg == -1)
    {
        char *poly_string = (char *)malloc(2 * sizeof(char));
        assert(poly_string);
        poly_string[0] = '0';
        poly_string[1] = '\0';
        return poly_string;
    }
    int size_p = size_int_str(p);
    int size_d = size_int_str(f->deg);
    char *poly_string = (char *)malloc(((size_p + size_d + 6) * (f->deg + 1) + 1) * sizeof(char));
    assert(poly_string);
    int index = 0;
    for (int i = 0; i <= f->deg; i++)
    {
        if (index != 0)
        {
            poly_string[index + 0] = ' ';
            poly_string[index + 1] = '+';
            poly_string[index + 2] = ' ';
            index += 3;
        }
        int size_coeff = size_int_str(f->coeffs[i]);
        for (int j = 0; j < size_p - size_coeff; j++)
            poly_string[index + j] = ' ';
        index += size_p - size_coeff;
        str_add_int(poly_string, &index, f->coeffs[i]);
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

void poly_print(const poly f)
{
    char *poly_string = str_poly_iso_length(f);
    printf("%s\n", poly_string);
    free(poly_string);
}

/******************************************************/

poly poly_new(void)
{
    poly f = (poly)malloc(sizeof(poly));
    assert(f);
    poly_set(f, -1, NULL);
    return f;
}

poly poly_new_set(int deg, ...)
{
    poly f = poly_new();
    if (deg > -1)
    {
        array coeffs = array_new(deg + 1);
        va_list valist;
        va_start(valist, deg);
        for (int i = 0; i <= deg; i++)
            coeffs[i] = va_arg(valist, int);
        va_end(valist);
        poly_set(f, deg, coeffs);
    }
    return f;
}

int next_coeff(char *str, array i)
{
    while (str[*i] != '\0' && str[*i] != '-' && (str[*i] < '0' || str[*i] > '9') && str[*i] != 'x' && str[*i] != 'X')
        *i += 1;
    if (str[*i] == '\0')
    {
        *i = -10;
        return -1;
    }
    int sign = 1;
    if (str[*i] == '-')
    {
        sign = -1;
        *i += 1;
    }
    while (str[*i] == ' ')
        *i += 1;
    if (str[*i] < '0' || str[*i] > '9')
        return sign;
    int n = 0;
    for (; str[*i] >= '0' && str[*i] <= '9'; *i += 1)
        n = 10 * n + str[*i] - '0';
    return sign * n;
}

array coeffs_from_str(char *str, array deg)
{
    int i = 0;
    int t = 100;
    array coeffs = array_new_zeros(t);
    int exp;
    for (int coeff = next_coeff(str, &i); i != -10; coeff = next_coeff(str, &i))
    {
        while (str[i] != '\0' && str[i] != 'x' && str[i] != 'X' && str[i] != '+' && str[i] != '-')
            i++;
        if (str[i] == '\0' || str[i] == '+' || str[i] == '-')
        {
            coeffs[0] += coeff;
            coeffs[0] = ((coeffs[0] % p) + p) % p;
            if (*deg < 0 && coeffs[0] != 0)
                *deg = 0;
            if (*deg == 0 && coeffs[0] == 0)
                *deg = -1;
        }
        else
        {
            while (str[i] != '\0' && str[i] != '*' && str[i] != '^' && str[i] != '+' && str[i] != '-')
                i++;
            if (str[i] == '\0' || str[i] == '+' || str[i] == '-')
                exp = 1;
            else
                exp = next_coeff(str, &i);
            if (i != -10 && exp >= 0)
            {
                if (t <= exp)
                {
                    array new_coeffs = array_new_zeros(2 * exp);
                    for (int j = 0; j < t; j++)
                        new_coeffs[j] = coeffs[j];
                    array_free(coeffs);
                    coeffs = new_coeffs;
                    t = 2 * exp;
                }
                coeffs[exp] += coeff;
                coeffs[exp] = ((coeffs[exp] % p) + p) % p;
                if (*deg < exp && coeffs[exp] != 0)
                    *deg = exp;
                if (*deg == exp && coeffs[exp] == 0)
                    while (*deg >= 0 && coeffs[*deg] == 0)
                        (*deg)--;
            }
        }
    }
    return coeffs;
}

poly poly_new_str(char *str)
{
    int deg = -1;
    array coeffs_str = coeffs_from_str(str, &deg);
    if (deg == -1)
    {
        array_free(coeffs_str);
        return poly_new();
    }
    poly f = poly_new();
    array coeffs = array_new(deg + 1);
    for (int i = 0; i <= deg; i++)
        coeffs[i] = coeffs_str[i];
    array_free(coeffs_str);
    poly_set(f, deg, coeffs);
    return f;
}

/******************************************************/

void poly_clear(poly f)
{
    array_free(f->coeffs);
    poly_set(f, -1, NULL);
}

void poly_free(poly f)
{
    poly_clear(f);
    free(f);
}

/******************************************************/

bool poly_equal(const poly f, const poly g)
{
    if (poly_deg(f) != poly_deg(g))
        return false;
    for (int i = 0; i < poly_deg(f); i++)
        if (poly_coeffs(f)[i] != poly_coeffs(g)[i])
            return false;
    return true;
}

poly poly_copy(const poly src)
{
    poly copy = poly_new();
    int deg = src->deg;
    if (deg > -1)
    {
        array coeffs = array_new(deg + 1);
        memcpy(coeffs, src->coeffs, (deg + 1) * sizeof(int));
        poly_set(copy, deg, coeffs);
    }
    return copy;
}

poly poly_rand(int deg)
{
    poly f = poly_new();
    array coeffs = array_new(deg + 1);
    for (int i = 0; i < deg; i++)
        coeffs[i] = zp_rand();
    while (coeffs[deg] == 0)
        coeffs[deg] = zp_rand();
    poly_set(f, deg, coeffs);
    return f;
}

void poly_rev(poly rop, const poly f)
{
    int deg_rev, *coeffs_rev;
    deg_rev = f->deg;
    for (int i = 0; f->coeffs[i] == 0; i++)
        deg_rev--;
    coeffs_rev = array_new(deg_rev + 1);
    for (int i = 0; i <= deg_rev; i++)
        coeffs_rev[i] = f->coeffs[f->deg - i];
    poly_clear(rop);
    poly_set(rop, deg_rev, coeffs_rev);
}

void poly_deriv(poly rop, const poly op)
{
    int deg, *coeffs;
    if (op->deg < 1)
    {
        deg = -1;
        coeffs = NULL;
    }
    else
    {
        deg = op->deg - 1;
        coeffs = array_new(deg + 1);
        for (int i = 0; i <= deg; i++)
            coeffs[i] = zp_mul(i + 1, op->coeffs[i + 1]);
    }
    poly_clear(rop);
    poly_set(rop, deg, coeffs);
}

/******************************************************/

void poly_add(poly rop, const poly op1, const poly op2)
{
    int deg, *coeffs;
    if (op1->deg > op2->deg)
    {
        deg = op1->deg;
        coeffs = array_new(deg + 1);
        int i;
        for (i = 0; i <= op2->deg; i++)
            coeffs[i] = zp_add(op1->coeffs[i], op2->coeffs[i]);
        for (; i <= deg; i++)
            coeffs[i] = op1->coeffs[i];
    }
    else if (op1->deg < op2->deg)
    {
        deg = op2->deg;
        coeffs = array_new(deg + 1);
        int i;
        for (i = 0; i <= op1->deg; i++)
            coeffs[i] = zp_add(op1->coeffs[i], op2->coeffs[i]);
        for (; i <= deg; i++)
            coeffs[i] = op2->coeffs[i];
    }
    else if (op1->deg == -1)
    {
        deg = -1;
        coeffs = NULL;
    }
    else
    {
        int i;
        for (i = op1->deg; i >= 0 && zp_add(op1->coeffs[i], op2->coeffs[i]) == 0; i--)
            ;
        if (i == -1)
        {
            deg = -1;
            coeffs = NULL;
        }
        else
        {
            deg = i;
            coeffs = array_new(deg + 1);
            for (int j = 0; j <= deg; j++)
                coeffs[j] = zp_add(op1->coeffs[j], op2->coeffs[j]);
        }
    }
    poly_clear(rop);
    poly_set(rop, deg, coeffs);
}

void poly_sub(poly rop, const poly op1, const poly op2)
{
    int deg, *coeffs;
    if (op1->deg > op2->deg)
    {
        deg = op1->deg;
        coeffs = array_new(deg + 1);
        int i;
        for (i = 0; i <= op2->deg; i++)
            coeffs[i] = zp_sub(op1->coeffs[i], op2->coeffs[i]);
        for (; i <= deg; i++)
            coeffs[i] = op1->coeffs[i];
    }
    else if (op1->deg < op2->deg)
    {
        deg = op2->deg;
        coeffs = array_new(deg + 1);
        int i;
        for (i = 0; i <= op1->deg; i++)
            coeffs[i] = zp_sub(op1->coeffs[i], op2->coeffs[i]);
        for (; i <= deg; i++)
            coeffs[i] = zp_opp(op2->coeffs[i]);
    }
    else if (op1->deg == -1)
    {
        deg = -1;
        coeffs = NULL;
    }
    else
    {
        int i;
        for (i = op1->deg; i >= 0 && zp_sub(op1->coeffs[i], op2->coeffs[i]) == 0; i--)
            ;
        if (i == -1)
        {
            deg = -1;
            coeffs = NULL;
        }
        else
        {
            deg = i;
            coeffs = array_new(deg + 1);
            for (int j = 0; j <= deg; j++)
                coeffs[j] = zp_sub(op1->coeffs[j], op2->coeffs[j]);
        }
    }
    poly_clear(rop);
    poly_set(rop, deg, coeffs);
}

void poly_mul(poly rop, const poly op1, const poly op2)
{
    int deg, *coeffs;
    if (op1->deg == -1 || op2->deg == -1)
    {
        deg = -1;
        coeffs = NULL;
    }
    else
    {
        deg = op1->deg + op2->deg;
        coeffs = array_new_zeros(deg + 1);
        for (int i = 0; i <= op1->deg; i++)
            for (int j = 0; j <= op2->deg; j++)
                coeffs[i + j] = zp_add(coeffs[i + j], zp_mul(op1->coeffs[i], op2->coeffs[j]));
    }
    poly_clear(rop);
    poly_set(rop, deg, coeffs);
}

void poly_mul_scalar(poly rop, int op1, const poly op2)
{
    int deg, *coeffs;
    if (op1 == 0 || op2->deg == -1)
    {
        deg = -1;
        coeffs = NULL;
    }
    else
    {
        deg = op2->deg;
        coeffs = array_new(deg + 1);
        for (int i = 0; i <= deg; i++)
            coeffs[i] = zp_mul(op2->coeffs[i], op1);
    }
    poly_clear(op2);
    poly_set(rop, deg, coeffs);
}

void poly_euc_div(poly q, poly r, const poly op1, const poly op2)
{
    int deg_q, *coeffs_q;
    int deg_r, *coeffs_r;
    array rest = array_new(op1->deg + 1);
    for (int i = 0; i <= op1->deg; i++)
        rest[i] = op1->coeffs[i];
    deg_q = op1->deg - op2->deg;
    coeffs_q = array_new_zeros(deg_q + 1);
    deg_r = op1->deg;
    int inv_LC_op2 = zp_inv(op2->coeffs[op2->deg]);
    while (deg_r >= op2->deg)
    {
        int c = zp_mul(rest[deg_r], inv_LC_op2);
        coeffs_q[deg_r - op2->deg] = c;
        for (int i = 0; i <= op2->deg; i++)
            rest[i + deg_r - op2->deg] = zp_sub(rest[i + deg_r - op2->deg], zp_mul(c, op2->coeffs[i]));
        while (deg_r >= 0 && rest[deg_r] == 0)
            deg_r--;
    }
    if (deg_r == -1)
        coeffs_r = NULL;
    else
        coeffs_r = array_new(deg_r + 1);
    for (int i = 0; i <= deg_r; i++)
        coeffs_r[i] = rest[i];
    array_free(rest);
    poly_set(q, deg_q, coeffs_q);
    poly_set(r, deg_r, coeffs_r);
}

void poly_xgcd(poly *d, poly *u, poly *v, const poly op1, const poly op2)
{
    *d = poly_copy(op1);
    *u = poly_new_str("1");
    *v = poly_new();
    poly r1 = poly_copy(op2);
    poly u1 = poly_new();
    poly v1 = poly_new_str("1");
    while (r1->deg > -1)
    {
        poly q = poly_new();
        poly r2 = poly_new();
        poly_euc_div(q, r2, *d, r1);
        poly u2 = poly_new();
        poly_mul(u2, q, u1);
        poly_sub(u2, *u, u2);
        poly v2 = poly_new();
        poly_mul(v2, q, v1);
        poly_sub(v2, *v, v2);
        poly_free(q);
        poly_free(*d);
        poly_free(*u);
        poly_free(*v);
        *d = r1;
        *u = u1;
        *v = v1;
        r1 = r2;
        u1 = u2;
        v1 = v2;
    }
    // We want d to be monic.
    int inv_LC_d = zp_inv((*d)->coeffs[(*d)->deg]);
    poly_mul_scalar(*d, inv_LC_d, *d);
    poly_mul_scalar(*u, inv_LC_d, *u);
    poly_mul_scalar(*v, inv_LC_d, *v);
    poly_free(r1);
    poly_free(u1);
    poly_free(v1);
}

void poly_xgcd_partial(poly *d, poly *u, poly *v, const poly op1, const poly op2, int limit)
{
    *d = poly_copy(op1);
    *u = poly_new_set(0, 1);
    *v = poly_new();
    poly r1 = poly_copy(op2);
    poly u1 = poly_new();
    poly v1 = poly_new_set(0, 1);
    while (limit <= (*d)->deg)
    {
        poly q = poly_new();
        poly r2 = poly_new();
        poly_euc_div(q, r2, *d, r1);
        poly u2 = poly_new();
        poly_mul(u2, q, u1);
        poly_sub(u2, *u, u2);
        poly v2 = poly_new();
        poly_mul(v2, q, v1);
        poly_sub(v2, *v, v2);
        poly_free(q);
        poly_free(*d);
        poly_free(*u);
        poly_free(*v);
        *d = r1;
        *u = u1;
        *v = v1;
        r1 = r2;
        u1 = u2;
        v1 = v2;
    }
    // We want d to be monic.
    poly_free(r1);
    poly_free(u1);
    poly_free(v1);
}

/******************************************************/

int poly_eval(const poly f, int x)
{
    int y = 0;
    int x_i = 1;
    for (int i = 0; i <= f->deg; i++)
    {
        y = zp_add(y, zp_mul(x_i, f->coeffs[i]));
        x_i = zp_mul(x_i, x);
    }
    return y;
}

array poly_eval_array(const poly f, const array tab, int tab_size)
{
    array eval = array_new(tab_size);
    for (int i = 0; i < tab_size; i++)
        eval[i] = poly_eval(f, tab[i]);
    return eval;
}

array poly_dft(const poly f)
{
    array eval = array_new(n);
    for (int i = 0; i < n; i++)
    {
        eval[i] = poly_eval(f, omegas[i]);
    }
    return eval;
}

void interpolation(poly rop, const array points, const array eval, int size)
{
    poly_clear(rop);
    poly f = poly_new_set(0, 1);
    for (int i = 0; i < size; i++)
    {
        poly x_m_ai = poly_new_set(1, zp_opp(points[i]), 1);
        poly_mul(f, f, x_m_ai);
        poly_free(x_m_ai);
    }
    poly fd = poly_new();
    poly_deriv(fd, f);
    for (int i = 0; i < size; i++)
    {
        poly x_m_ai = poly_new_set(1, zp_opp(points[i]), 1);
        int fd_ai = poly_eval(fd, points[i]);
        poly_mul_scalar(x_m_ai, fd_ai, x_m_ai);
        poly l_i = poly_new();
        poly r = poly_new();
        poly_euc_div(l_i, r, f, x_m_ai);
        poly_mul_scalar(l_i, eval[i], l_i);
        poly_add(rop, rop, l_i);
        poly_free(x_m_ai);
        poly_free(l_i);
        poly_free(r);
    }
    poly_free(f);
    poly_free(fd);
}

void poly_inv_dft(poly rop, array eval)
{
    interpolation(rop, omegas, eval, n);
}

/******************************************************/

array poly_fft(const poly f)
{
    // We need some padding.
    array coeffs = array_new_zeros(n);
    for (int i = 0; i <= f->deg; i++)
        coeffs[i] = f->coeffs[i];
    array eval = array_fft(coeffs, n);
    array_free(coeffs);
    return eval;
}

void poly_inv_fft(poly rop, array eval)
{
    poly_clear(rop);
    array f_coeffs = array_inv_fft(eval, n);
    int deg = n - 1;
    while (deg >= 0 && f_coeffs[deg] == 0)
        deg--;
    array coeffs;
    if (deg == -1)
        coeffs = NULL;
    else
    {
        coeffs = array_new(deg + 1);
        for (int i = 0; i <= deg; i++)
            coeffs[i] = f_coeffs[i];
    }
    array_free(f_coeffs);
    poly_set(rop, deg, coeffs);
}

/******************************************************/

void poly_fast_mul(poly rop, const poly op1, const poly op2)
{
    poly_clear(rop);
    array eval_op1 = poly_fft(op1);
    array eval_op2 = poly_fft(op2);
    array eval_rop = array_new(n);
    for (int i = 0; i < n; i++)
        eval_rop[i] = zp_mul(eval_op1[i], eval_op2[i]);
    poly_inv_fft(rop, eval_rop);
    array_free(eval_op1);
    array_free(eval_op2);
    array_free(eval_rop);
}

void poly_fast_euc_div(poly quo, poly rem, const poly op1, const poly op2)
{
    int deg_p = op1->deg;
    int deg_d = op2->deg;
    array coeffs_quo = NULL;
    array coeffs_rem = NULL;
    formal_serie_fast_euc_div(&coeffs_quo, &coeffs_rem, op1->coeffs, deg_p, op2->coeffs, deg_d);
    poly_set(quo, deg_p - deg_d, coeffs_quo);
    int deg_rem = deg_d - 1;
    while (deg_rem >= 0 && coeffs_rem[deg_rem] == 0)
        deg_rem--;
    if (deg_rem == -1)
        poly_set(rem, -1, NULL);
    else
    {
        array coeffs = array_new(deg_rem + 1);
        for (int i = 0; i <= deg_rem; i++)
            coeffs[i] = coeffs_rem[i];
        poly_set(rem, deg_rem, coeffs);
    }
    array_free(coeffs_rem);
}
