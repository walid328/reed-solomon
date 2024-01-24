#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "assert.h"
#include "polynomial.h"
#include "finite_field.h"
#include "field_spec.h"

struct polynomial
{
    int degree;
    int *coeffs;
};

/******************************************************/

void poly_set(poly *f, int deg, int *coeffs)
{
    f->degree = deg;
    f->coeffs = coeffs;
}

int poly_degree(poly *f)
{
    return f->degree;
}

int *poly_coeffs(poly *f)
{
    return f->coeffs;
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
        ASSERT(poly_string);
        poly_string[0] = '0';
        poly_string[1] = '\0';
        return poly_string;
    }
    int size_p = size_int_str(p);
    int size_d = size_int_str(q->degree);
    char *poly_string = (char *)malloc(((size_p + size_d + 6) * (q->degree + 1) + 1) * sizeof(char));
    ASSERT(poly_string);
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

void poly_print(poly *q)
{
    char *poly_string = str_poly_iso_length(q);
    printf("%s\n", poly_string);
    free(poly_string);
}

/******************************************************/

poly *poly_new(void)
{
    poly *f = (poly *)malloc(sizeof(poly));
    ASSERT(f);
    poly_set(f, -1, NULL);
    return f;
}

poly *poly_new_1(void)
{
    poly *f = poly_new();
    int *coeffs = (int *)malloc(1 * sizeof(int));
    ASSERT(coeffs);
    coeffs[0] = 1;
    poly_set(f, 0, coeffs);
    return f;
}

poly *poly_new_from_coeffs(int deg, ...)
{
    if (deg == -1)
        return poly_new();
    poly *f = poly_new();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    ASSERT(coeffs);
    va_list valist;
    va_start(valist, deg);
    for (int i = 0; i <= deg; i++)
        coeffs[i] = va_arg(valist, int);
    va_end(valist);
    poly_set(f, deg, coeffs);
    return f;
}

int next_coeff(char *str, int *i)
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

int *coeffs_from_str(char *str, int *deg)
{
    int i = 0;
    int t = 100;
    int *coeffs = (int *)calloc(t, sizeof(int));
    ASSERT(coeffs);
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
                    int *new_coeffs = (int *)calloc(2 * exp, sizeof(int));
                    ASSERT(new_coeffs);
                    for (int j = 0; j < t; j++)
                        new_coeffs[j] = coeffs[j];
                    free(coeffs);
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

poly *poly_new_from_str(char *str)
{
    int deg = -1;
    int *coeffs_str = coeffs_from_str(str, &deg);
    if (deg == -1)
    {
        free(coeffs_str);
        return poly_new();
    }
    poly *f = poly_new();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    ASSERT(coeffs);
    for (int i = 0; i <= deg; i++)
        coeffs[i] = coeffs_str[i];
    free(coeffs_str);
    poly_set(f, deg, coeffs);
    return f;
}

poly *poly_new_from_copy(poly *source)
{
    int deg = source->degree;
    if (deg == -1)
        return poly_new();
    poly *copy = poly_new();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    ASSERT(coeffs);
    for (int i = 0; i <= deg; i++)
        coeffs[i] = source->coeffs[i];
    poly_set(copy, deg, coeffs);
    return copy;
}

poly *poly_new_rand(int deg)
{
    poly *f = poly_new();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    ASSERT(coeffs);
    for (int i = 0; i <= deg; i++)
        coeffs[i] = zp_rand();
    while (coeffs[deg] == 0)
        coeffs[deg] = zp_rand();
    poly_set(f, deg, coeffs);
    return f;
}

/******************************************************/

void poly_clear(poly *f)
{
    f->degree = -1;
    f->coeffs = NULL;
}

void poly_clear_full(poly *f)
{
    f->degree = -1;
    free(f->coeffs);
    f->coeffs = NULL;
}

void poly_free(poly *f)
{
    free(f);
}

void poly_free_full(poly *f)
{
    free(f->coeffs);
    free(f);
}

/******************************************************/

void poly_add(poly *rop, poly *op1, poly *op2)
{
    int deg;
    int *coeffs;
    if (op1->degree != op2->degree)
    {
        poly *max = op1->degree > op2->degree ? op1 : op2;
        poly *min = op1->degree > op2->degree ? op2 : op1;
        deg = max->degree;
        coeffs = (int *)malloc((deg + 1) * sizeof(poly));
        ASSERT(coeffs);
        int i;
        for (i = 0; i <= min->degree; i++)
            coeffs[i] = zp_add(op1->coeffs[i], op2->coeffs[i]);
        for (; i <= deg; i++)
            coeffs[i] = max->coeffs[i];
    }
    else if (op1->degree == -1)
    {
        deg = -1;
        coeffs = NULL;
    }
    else
    {
        int i;
        for (i = op1->degree; i >= 0 && zp_add(op1->coeffs[i], op2->coeffs[i]) == 0; i--)
            ;
        if (i == -1)
        {
            deg = -1;
            coeffs = NULL;
        }
        else
        {
            deg = i;
            coeffs = (int *)malloc((deg + 1) * sizeof(int));
            ASSERT(coeffs);
            for (int j = 0; j <= deg; j++)
                coeffs[j] = zp_add(op1->coeffs[j], op2->coeffs[j]);
        }
    }
    if (rop == op1 || rop == op2)
        poly_clear_full(rop);
    poly_set(rop, deg, coeffs);
}

void poly_sub(poly *rop, poly *op1, poly *op2)
{
    int deg;
    int *coeffs;
    if (op1->degree > op2->degree)
    {
        deg = op1->degree;
        coeffs = (int *)malloc((deg + 1) * sizeof(int));
        ASSERT(coeffs);
        int i;
        for (i = 0; i <= op2->degree; i++)
            coeffs[i] = zp_sub(op1->coeffs[i], op2->coeffs[i]);
        for (; i <= deg; i++)
            coeffs[i] = op1->coeffs[i];
    }
    else if (op1->degree < op2->degree)
    {
        deg = op2->degree;
        coeffs = (int *)malloc((deg + 1) * sizeof(int));
        ASSERT(coeffs);
        int i;
        for (i = 0; i <= op1->degree; i++)
            coeffs[i] = zp_sub(op1->coeffs[i], op2->coeffs[i]);
        for (; i <= deg; i++)
            coeffs[i] = zp_opp(op2->coeffs[i]);
    }
    else if (op1->degree == -1)
    {
        deg = -1;
        coeffs = NULL;
    }
    else
    {
        int i;
        for (i = op1->degree; i >= 0 && zp_sub(op1->coeffs[i], op2->coeffs[i]) == 0; i--)
            ;
        if (i == -1)
        {
            deg = -1;
            coeffs = NULL;
        }
        else
        {
            deg = i;
            coeffs = (int *)malloc((deg + 1) * sizeof(int));
            ASSERT(coeffs);
            for (int j = 0; j <= deg; j++)
                coeffs[j] = zp_sub(op1->coeffs[j], op2->coeffs[j]);
        }
    }
    if (rop == op1 || rop == op2)
        poly_clear_full(rop);
    poly_set(rop, deg, coeffs);
}

void poly_mul(poly *rop, poly *op1, poly *op2)
{
    int deg;
    int *coeffs;
    if (op1->degree == -1 || op2->degree == -1)
    {
        deg = -1;
        coeffs = NULL;
    }
    else
    {
        deg = op1->degree + op2->degree;
        coeffs = (int *)calloc(deg + 1, sizeof(int));
        ASSERT(coeffs);
        for (int i = 0; i <= op1->degree; i++)
            for (int j = 0; j <= op2->degree; j++)
                coeffs[i + j] = zp_add(coeffs[i + j], zp_mul(op1->coeffs[i], op2->coeffs[j]));
    }
    if (rop == op1 || rop == op2)
        poly_clear_full(rop);
    poly_set(rop, deg, coeffs);
}

void poly_mul_scalar(poly *rop, int op1, poly *op2)
{
    int deg;
    int *coeffs;
    if (op1 == 0 || op2->degree == -1)
    {
        deg = -1;
        coeffs = NULL;
    }
    else
    {
        deg = op2->degree;
        coeffs = (int *)malloc((deg + 1) * sizeof(int));
        ASSERT(coeffs);
        for (int i = 0; i <= deg; i++)
            coeffs[i] = zp_mul(op2->coeffs[i], op1);
    }
    if (rop == op2)
        poly_clear_full(op2);
    poly_set(rop, deg, coeffs);
}

void poly_euc_div(poly *q, poly *r, poly *op1, poly *op2)
{
    int deg_q, deg_r;
    int *coeffs_q, *coeffs_r;
    int *rest = (int *)malloc((op1->degree + 1) * sizeof(int));
    for (int i = 0; i <= op1->degree; i++)
        rest[i] = op1->coeffs[i];
    deg_q = op1->degree - op2->degree;
    coeffs_q = (int *)calloc(deg_q + 1, sizeof(int));
    ASSERT(coeffs_q);
    deg_r = op1->degree;
    int inv_LC_op2 = zp_inv(op2->coeffs[op2->degree]);
    while (deg_r >= op2->degree)
    {
        int c = zp_mul(rest[deg_r], inv_LC_op2);
        coeffs_q[deg_r - op2->degree] = c;
        for (int i = 0; i <= op2->degree; i++)
            rest[i + deg_r - op2->degree] = zp_sub(rest[i + deg_r - op2->degree], zp_mul(c, op2->coeffs[i]));
        while (deg_r >= 0 && rest[deg_r] == 0)
            deg_r--;
    }
    if (deg_r == -1)
        coeffs_r = NULL;
    else
    {
        coeffs_r = (int *)malloc((deg_r + 1) * sizeof(int));
        ASSERT(coeffs_r);
    }
    for (int i = 0; i <= deg_r; i++)
        coeffs_r[i] = rest[i];
    free(rest);
    poly_set(q, deg_q, coeffs_q);
    poly_set(r, deg_r, coeffs_r);
}

void poly_xgcd(poly **d, poly **u, poly **v, poly *op1, poly *op2)
{
    *d = poly_new_from_copy(op1);
    *u = poly_new_1();
    *v = poly_new();
    poly *r1 = poly_new_from_copy(op2);
    poly *u1 = poly_new();
    poly *v1 = poly_new_1();
    while (r1->degree > -1)
    {
        poly *q = poly_new();
        poly *r2 = poly_new();
        poly_euc_div(q, r2, *d, r1);
        poly *u2 = poly_new();
        poly_mul(u2, q, u1);
        poly_sub(u2, *u, u2);
        poly *v2 = poly_new();
        poly_mul(v2, q, v1);
        poly_sub(v2, *v, v2);
        poly_free_full(q);
        poly_free_full(*d);
        poly_free_full(*u);
        poly_free_full(*v);
        *d = r1;
        *u = u1;
        *v = v1;
        r1 = r2;
        u1 = u2;
        v1 = v2;
    }
    // We want monic polynomials.
    int inv_LC_d = zp_inv((*d)->coeffs[(*d)->degree]);
    poly_mul_scalar(*d, inv_LC_d, *d);
    poly_mul_scalar(*u, inv_LC_d, *u);
    poly_mul_scalar(*v, inv_LC_d, *v);
    poly_free_full(r1);
    poly_free_full(u1);
    poly_free_full(v1);
}

void poly_xgcd_partial(poly **d, poly **u, poly **v, poly *op1, poly *op2, int n)
{
    *d = poly_new_from_copy(op1);
    *u = poly_new_1();
    *v = poly_new();
    poly *r1 = poly_new_from_copy(op2);
    poly *u1 = poly_new();
    poly *v1 = poly_new_1();
    while (n <= (*d)->degree)
    {
        poly *q = poly_new();
        poly *r2 = poly_new();
        poly_euc_div(q, r2, *d, r1);
        poly *u2 = poly_new();
        poly_mul(u2, q, u1);
        poly_sub(u2, *u, u2);
        poly *v2 = poly_new();
        poly_mul(v2, q, v1);
        poly_sub(v2, *v, v2);
        poly_free_full(q);
        poly_free_full(*d);
        poly_free_full(*u);
        poly_free_full(*v);
        *d = r1;
        *u = u1;
        *v = v1;
        r1 = r2;
        u1 = u2;
        v1 = v2;
    }
    poly_free_full(r1);
    poly_free_full(u1);
    poly_free_full(v1);
}

/******************************************************/

void poly_deriv(poly *rop, poly *op)
{
    int deg;
    int *coeffs;
    if (op->degree < 1)
    {
        deg = -1;
        coeffs = NULL;
    }
    else
    {
        deg = op->degree - 1;
        coeffs = (int *)malloc((deg + 1) * sizeof(int));
        ASSERT(coeffs);
        for (int i = 0; i <= deg; i++)
            coeffs[i] = zp_mul(i + 1, op->coeffs[i + 1]);
    }
    if (rop == op)
        poly_clear_full(op);
    poly_set(rop, deg, coeffs);
}

int poly_eval(poly *f, int x)
{
    int y = 0;
    int x_i = 1;
    for (int i = 0; i <= f->degree; i++)
    {
        y = zp_add(y, zp_mul(x_i, f->coeffs[i]));
        x_i = zp_mul(x_i, x);
    }
    return y;
}

void poly_eval_multi(poly *f, int n, int *a, int *b)
{
    for (int i = 0; i <= n; i++)
        b[i] = poly_eval(f, a[i]);
}

void poly_dft(poly *f, int **eval)
{
    int q, d;
    for (q = p - 1, d = 1; (q & 1) == 0; q >>= 1, d <<= 1)
        ;
    if (f->degree >= d)
    {
        fprintf(stderr, "poly_dft degree too high");
        exit(EXIT_FAILURE);
    }
    int omega = zp_prim_root_min(q, d);
    *eval = (int *)malloc(d * sizeof(int));
    ASSERT(*eval);
    int omega_i = 1;
    for (int i = 0; i < d; i++, omega_i = zp_mul(omega_i, omega))
        (*eval)[i] = poly_eval(f, omega_i);
}

poly *interpolation(int *a, int *b, int n)
{
    poly *f = poly_new_1();
    for (int i = 0; i < n; i++)
    {
        poly *x_m_ai = poly_new_from_coeffs(1, zp_opp(a[i]), 1);
        poly_mul(f, f, x_m_ai);
        poly_free_full(x_m_ai);
    }
    poly *fd = poly_new();
    poly_deriv(fd, f);
    poly *g = poly_new();
    for (int i = 0; i < n; i++)
    {
        poly *x_m_ai = poly_new_from_coeffs(1, zp_opp(a[i]), 1);
        int fd_ai = poly_eval(fd, a[i]);
        poly_mul_scalar(x_m_ai, fd_ai, x_m_ai);
        poly *l_i = poly_new();
        poly *r = poly_new();
        poly_euc_div(l_i, r, f, x_m_ai);
        poly_mul_scalar(l_i, b[i], l_i);
        poly_add(g, g, l_i);
        poly_free_full(x_m_ai);
        poly_free_full(l_i);
        poly_free_full(r);
    }
    poly_free_full(f);
    poly_free_full(fd);
    return g;
}

/******************************************************/

void array_split_all(int **even, int *even_size, int **odd, int *odd_size, int *tab, int tab_size)
{
    *even_size = tab_size / 2;
    *odd_size = *even_size;
    if (tab_size & 1)
        (*even_size)++;
    *even = (int *)malloc(*even_size * sizeof(int));
    ASSERT(*even);
    *odd = (int *)malloc(*odd_size * sizeof(int));
    ASSERT(*odd);
    int even_cnt = 0;
    int odd_cnt = 0;
    for (int i = 0; i < tab_size; i++)
        if (i & 1)
            (*odd)[odd_cnt++] = tab[i];
        else
            (*even)[even_cnt++] = tab[i];
}

void array_split(int **even, int **odd, int *tab, int tab_size)
{
    *even = (int *)malloc(tab_size / 2 * sizeof(int));
    ASSERT(even);
    *odd = (int *)malloc(tab_size / 2 * sizeof(int));
    ASSERT(odd);
    for (int i = 0; i < tab_size / 2; i++)
    {
        (*even)[i] = tab[2 * i];
        (*odd)[i] = tab[2 * i + 1];
    }
}

void array_merge(int **tab, int *even, int *odd, int subtab_size)
{
    *tab = (int *)malloc(2 * subtab_size * sizeof(int));
    ASSERT(*tab);
    for (int i = 0; i < subtab_size; i++)
    {
        (*tab)[2 * i] = even[i];
        (*tab)[2 * i + 1] = odd[i];
    }
}

void poly_split(poly *even, poly *odd, poly *f)
{
    int *coeffs_even = NULL;
    int even_size = 0;
    int *coeffs_odd = NULL;
    int odd_size = 0;
    array_split_all(&coeffs_even, &even_size, &coeffs_odd, &odd_size, f->coeffs, f->degree + 1);
    poly_set(even, even_size - 1, coeffs_even);
    poly_set(odd, odd_size - 1, coeffs_odd);
}

int *fft(int *f, int d)
{
    if (d == 1)
    {
        int *eval = (int *)malloc(sizeof(int));
        ASSERT(eval);
        eval[0] = f[0];
        return eval;
    }
    int *r0 = (int *)malloc((d / 2) * sizeof(int));
    ASSERT(r0);
    for (int i = 0; i < d / 2; i++)
        r0[i] = zp_add(f[i], f[i + d / 2]);
    int *r1 = (int *)malloc((d / 2) * sizeof(int));
    ASSERT(r1);
    for (int i = 0; i < d / 2; i++)
        r1[i] = zp_mul(zp_sub(f[i], f[i + d / 2]), omegas[i * (n / d)]);
    int *eval_r0 = fft(r0, d / 2);
    int *eval_r1 = fft(r1, d / 2);
    free(r0);
    free(r1);
    int *eval;
    array_merge(&eval, eval_r0, eval_r1, d / 2);
    free(eval_r0);
    free(eval_r1);
    return eval;
}

void poly_fft(poly *f, int **eval)
{
    int q, d;
    for (q = p - 1, d = 1; (q & 1) == 0; q >>= 1, d <<= 1)
        ;
    if (f->degree >= d)
    {
        fprintf(stderr, "poly_fft degree too high");
        exit(EXIT_FAILURE);
    }
    int *coeffs = (int *)calloc(d, sizeof(int));
    ASSERT(coeffs);
    for (int i = 0; i <= f->degree; i++)
        coeffs[i] = f->coeffs[i];
    // memcpy(coeffs, f->coeffs, (f->degree + 1) * sizeof(int));
    *eval = fft(coeffs, d);
    free(coeffs);
}

int *inv_fft(int *eval, int d)
{
    if (d == 1)
    {
        int *f = (int *)malloc(sizeof(int));
        ASSERT(f);
        f[0] = eval[0];
        return f;
    }
    int *eval_r0;
    int *eval_r1;
    array_split(&eval_r0, &eval_r1, eval, d);
    int *r0 = inv_fft(eval_r0, d / 2);
    int *r1 = inv_fft(eval_r1, d / 2);
    free(eval_r0);
    free(eval_r1);
    for (int i = 1; i <= d / 2; i++)
        r1[d / 2 - i] = zp_opp(zp_mul(r1[d / 2 - i], omegas[i * (n / d)]));
    int *f = (int *)malloc(d * sizeof(int));
    ASSERT(f);
    for (int i = 0; i < d / 2; i++)
    {
        f[i] = zp_mul(zp_add(r0[i], r1[i]), (p + 1) / 2);
        f[i + d / 2] = zp_mul(zp_sub(r0[i], r1[i]), (p + 1) / 2);
    }
    free(r0);
    free(r1);
    return f;
}

void poly_inv_fft(poly *f, int *eval)
{
    int *f_coeffs = inv_fft(eval, n);
    int degree = n - 1;
    while (degree >= 0 && f_coeffs[degree] == 0)
        degree--;
    int *coeffs;
    if (degree == -1)
    {
        coeffs = NULL;
    }
    else
    {
        coeffs = (int *)malloc((degree + 1) * sizeof(int));
        ASSERT(coeffs);
        for (int i = 0; i <= degree; i++)
            coeffs[i] = f_coeffs[i];
    }
    free(f_coeffs);
    poly_set(f, degree, coeffs);
}

int *formal_serie_mul(int *P, int *Q, int d)
{
    int *eval_P = fft(P, d);
    int *eval_Q = fft(Q, d);
    int *eval_res = (int *)malloc(d * sizeof(int));
    ASSERT(eval_res);
    for (int i = 0; i < d; i++)
        eval_res[i] = zp_mul(eval_P[i], eval_Q[i]);
    int *res = inv_fft(eval_res, d);
    free(eval_P);
    free(eval_Q);
    free(eval_res);
    return res;
}

int *formal_serie_inv(int *P, int d)
{
    if (d == 1)
    {
        int *Q = (int *)malloc(sizeof(int));
        ASSERT(Q);
        Q[0] = zp_inv(P[0]);
        return Q;
    }
    int *Q_0 = formal_serie_inv(P, d / 2);
    int *Q_0d = (int *)calloc(2 * d, sizeof(int));
    ASSERT(Q_0d);
    for (int i = 0; i < d / 2; i++)
        Q_0d[i] = Q_0[i];
    int *Q_a = formal_serie_mul(Q_0d, Q_0d, 2 * d);
    int *P_0 = (int *)calloc(2 * d, sizeof(int));
    ASSERT(P_0);
    for (int i = 0; i < d; i++)
        P_0[i] = P[i];
    int *Q_a2 = formal_serie_mul(P_0, Q_a, 2 * d);
    int *Q = (int *)calloc(2 * d, sizeof(int));
    ASSERT(Q);
    for (int i = 0; i < d / 2; i++)
        Q[i] = Q_0[i];
    for (int i = d / 2; i < d; i++)
        Q[i] = zp_opp(Q_a2[i]);
    free(Q_0);
    free(Q_0d);
    free(P_0);
    free(Q_a);
    free(Q_a2);
    return Q;
}

void poly_rev(poly *rop, poly *op)
{
    int degree_rop = op->degree;
    int *coeffs_rop = (int *)malloc(degree_rop * sizeof(int));
    ASSERT(coeffs_rop);
    for (int i = 0; i < degree_rop; i++)
        coeffs_rop[i] = op->coeffs[degree_rop - 1];
    poly_set(rop, degree_rop, coeffs_rop);
}

void poly_fast_mul(poly *rop, poly *op1, poly *op2)
{
    int *eval_op1 = NULL;
    poly_fft(op1, &eval_op1);
    int *eval_op2 = NULL;
    poly_fft(op2, &eval_op2);
    int *eval_rop = (int *)malloc(n * sizeof(int));
    ASSERT(eval_rop);
    for (int i = 0; i < n; i++)
        eval_rop[i] = zp_mul(eval_op1[i], eval_op2[i]);
    poly_inv_fft(rop, eval_rop);
    free(eval_op1);
    free(eval_op2);
    free(eval_rop);
}

void poly_fast_euc_div(poly *quo, poly *rem, poly *op1, poly *op2)
{
	int deg_p = op1->degree;
	int deg_d = op2->degree;
	int *quo_coeffs = NULL;
	int *res_coeffs	= NULL;
	formal_serie_fast_euc_div(&quo_coeffs, &res_coeffs, op1->coeffs, deg_p, op2->coeffs, deg_d);
	poly_set(quo, deg_p - deg_d, quo_coeffs);
	int degree_r = deg_d - 1;
	while (degree_r >= 0 && res_coeffs[degree_r] == 0)
		degree_r--;
	if (degree_r == -1)
		poly_set(rem, -1, NULL);
	else
	{
		int *r_coeffs = (int *)malloc((degree_r + 1) * sizeof(int));
		ASSERT(r_coeffs);
		for (int i = 0; i <= degree_r; i++)
			r_coeffs[i] = res_coeffs[i];
		poly_set(rem, degree_r, r_coeffs);
	}
	free(res_coeffs);
}

void formal_serie_fast_euc_div(int **quo, int **rem, int *op_p, int deg_p, int *op_d, int deg_d)
{
	// dm is the smallest power of 2 greater or equal to deg_d
	int dm = 1;
	while (dm < deg_d)
		dm <<= 1;
	int dnmm = 1;
	// dnmm is the smallest power of 2 greater or equal to op_p - deg_d + 1
	while (dnmm < deg_p - deg_d + 1)
		dnmm <<= 1;
	int *op_pp = (int *)calloc(2 * dnmm, sizeof(int));
	ASSERT(op_pp);
	int *op_dp = (int *)calloc(dnmm, sizeof(int));
	ASSERT(op_dp);
	for (int i = 0; i < deg_p - deg_d + 1; i++)
		op_pp[i] = op_p[deg_p - i];
	for (int i = 0; (i < (deg_p - deg_d + 1)) && (i < (deg_d + 1)); i++)
		op_dp[i] = op_d[deg_d - i];
	int *op_dpp = formal_serie_inv(op_dp, dnmm);
	int *qp = formal_serie_mul(op_pp, op_dpp, 2 * dnmm);
	*quo = (int *)malloc(deg_p - deg_d + 1 * sizeof(int));
	ASSERT(*quo);
	for (int i = 0; i < deg_p - deg_d + 1; i++)
		(*quo)[i] = qp[deg_p - deg_d - i];
	int *qdm = (int *)calloc(2 * dm, sizeof(int));
	ASSERT(qdm);
	for (int i = 0; i < deg_p - deg_d + 1 && i < deg_d; i++)
		qdm[i] = (*quo)[i];
	int *op_ddm = (int *)calloc(2 * dm, sizeof(int));
	ASSERT(op_ddm);
	for (int i = 0; i < deg_d; i++)
		op_ddm[i] = op_d[i];
	int *qd = formal_serie_mul(qdm, op_ddm, 2 * dm);
	*rem = (int *)malloc(deg_d * sizeof(int));
	ASSERT(*rem);
	for (int i = 0; i < deg_d; i++)
		(*rem)[i] = zp_sub(op_p[i], qd[i]);
	free(op_pp);
	free(op_dp);
	free(op_dpp);
	free(qp);
	free(qdm);
	free(op_ddm);
	free(qd);
}