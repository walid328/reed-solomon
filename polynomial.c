#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>

#include "polynomial.h"
#include "finite_field.h"

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

int degree(poly *f)
{
    return f->degree;
}

int *coeffs(poly *f)
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
        assert(poly_string);
        poly_string[0] = '0';
        poly_string[1] = '\0';
        return poly_string;
    }
    int size_p = size_int_str(p);
    int size_d = size_int_str(q->degree);
    char *poly_string = (char *)malloc(((size_p + size_d + 6) * (q->degree + 1) + 1) * sizeof(char));
    assert(poly_string);
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
    assert(f);
    poly_set(f, -1, NULL);
    return f;
}

poly *poly_new_0(void)
{
    poly *f = poly_new();
    int *coeffs = (int *)calloc(1, sizeof(int));
    assert(coeffs);
    poly_set(f, -1, coeffs);
    return f;
}

poly *poly_new_1(void)
{
    poly *f = poly_new();
    int *coeffs = (int *)malloc(1 * sizeof(int));
    assert(coeffs);
    coeffs[0] = 1;
    poly_set(f, 0, coeffs);
    return f;
}

poly *poly_new_from_coeffs(int deg, ...)
{
    if (deg == -1)
        return poly_new_0();
    poly *f = poly_new();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    assert(coeffs);
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
    assert(coeffs);
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
                    assert(new_coeffs);
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
        return poly_new_0();
    }
    poly *f = poly_new();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    assert(coeffs);
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
        return poly_new_0();
    poly *copy = poly_new();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    assert(coeffs);
    for (int i = 0; i <= deg; i++)
        coeffs[i] = source->coeffs[i];
    poly_set(copy, deg, coeffs);
    return copy;
}

poly *poly_new_rand(int deg)
{
    poly *f = poly_new();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    assert(coeffs);
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
    assert(f->coeffs);
    free(f->coeffs);
    f->coeffs = NULL;
}

void poly_free(poly *f)
{
    assert(f);
    free(f);
}

void poly_free_full(poly *f)
{
    assert(f);
    assert(f->coeffs);
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
        assert(coeffs);
        int i;
        for (i = 0; i <= min->degree; i++)
            coeffs[i] = zp_add(op1->coeffs[i], op2->coeffs[i]);
        for (; i <= deg; i++)
            coeffs[i] = max->coeffs[i];
    }
    else if (op1->degree == -1)
    {
        deg = -1;
        coeffs = (int *)calloc(1, sizeof(int));
        assert(coeffs);
    }
    else
    {
        int i;
        for (i = op1->degree; i >= 0 && zp_add(op1->coeffs[i], op2->coeffs[i]) == 0; i--)
            ;
        if (i == -1)
        {
            deg = -1;
            coeffs = (int *)calloc(1, sizeof(int));
            assert(coeffs);
        }
        else
        {
            deg = i;
            coeffs = (int *)malloc((deg + 1) * sizeof(int));
            assert(coeffs);
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
        assert(coeffs);
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
        assert(coeffs);
        int i;
        for (i = 0; i <= op1->degree; i++)
            coeffs[i] = zp_sub(op1->coeffs[i], op2->coeffs[i]);
        for (; i <= deg; i++)
            coeffs[i] = zp_opp(op2->coeffs[i]);
    }
    else if (op1->degree == -1)
    {
        deg = -1;
        coeffs = (int *)calloc(1, sizeof(int));
        assert(coeffs);
    }
    else
    {
        int i;
        for (i = op1->degree; i >= 0 && zp_sub(op1->coeffs[i], op2->coeffs[i]) == 0; i--)
            ;
        if (i == -1)
        {
            deg = -1;
            coeffs = (int *)calloc(1, sizeof(int));
            assert(coeffs);
        }
        else
        {
            deg = i;
            coeffs = (int *)malloc((deg + 1) * sizeof(int));
            assert(coeffs);
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
        coeffs = (int *)calloc(1, sizeof(int));
        assert(coeffs);
    }
    else
    {
        deg = op1->degree + op2->degree;
        coeffs = (int *)calloc(deg + 1, sizeof(int));
        assert(coeffs);
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
        coeffs = (int *)calloc(1, sizeof(int));
        assert(coeffs);
    }
    else
    {
        deg = op2->degree;
        coeffs = (int *)malloc((deg + 1) * sizeof(int));
        assert(coeffs);
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
    assert(coeffs_q);
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
        coeffs_r = (int *)calloc(1, sizeof(int));
    else
        coeffs_r = (int *)malloc((deg_r + 1) * sizeof(int));
    assert(coeffs_r);
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
    *v = poly_new_0();
    poly *r1 = poly_new_from_copy(op2);
    poly *u1 = poly_new_0();
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
    *v = poly_new_0();
    poly *r1 = poly_new_from_copy(op2);
    poly *u1 = poly_new_0();
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
        coeffs = (int *)calloc(1, sizeof(int));
        assert(coeffs);
    }
    else
    {
        deg = op->degree - 1;
        coeffs = (int *)malloc((deg + 1) * sizeof(int));
        assert(coeffs);
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
    assert(*eval);
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
    poly *g = poly_new_0();
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
    assert(*even);
    *odd = (int *)malloc(*odd_size * sizeof(int));
    assert(*odd);
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
    assert(even);
    *odd = (int *)malloc(tab_size / 2 * sizeof(int));
    assert(odd);
    for (int i = 0; i < tab_size / 2; i++)
    {
        (*even)[i] = tab[2 * i];
        (*odd)[i] = tab[2 * i + 1];
    }
}

void array_merge(int **tab, int *even, int *odd, int subtab_size)
{
    *tab = (int *)malloc(2 * subtab_size * sizeof(int));
    assert(*tab);
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

int *fft(int *f, int d, int omega)
{
    if (d == 1)
    {
        int *eval = (int *)malloc(sizeof(int));
        assert(eval);
        eval[0] = f[0];
        return eval;
    }
    int *r0 = (int *)malloc((d / 2) * sizeof(int));
    assert(r0);
    for (int i = 0; i < d / 2; i++)
        r0[i] = zp_add(f[i], f[i + d / 2]);
    int *r1 = (int *)malloc((d / 2) * sizeof(int));
    assert(r1);
    int omega_i = 1;
    for (int i = 0; i < d / 2; i++)
    {
        r1[i] = zp_mul(zp_sub(f[i], f[i + d / 2]), omega_i);
        omega_i = zp_mul(omega_i, omega);
    }
    int *eval_r0 = fft(r0, d / 2, zp_mul(omega, omega));
    int *eval_r1 = fft(r1, d / 2, zp_mul(omega, omega));
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
    int omega = zp_prim_root_min(q, d);
    int *coeffs = (int *)calloc(d, sizeof(int));
    assert(coeffs);
    for (int i = 0; i <= f->degree; i++)
        coeffs[i] = f->coeffs[i];
    // memcpy(coeffs, f->coeffs, (f->degree + 1) * sizeof(int));
    *eval = fft(coeffs, d, omega);
    free(coeffs);
}

int *inv_fft(int *eval, int d, int omega)
{
    if (d == 1)
    {
        int *f = (int *)malloc(sizeof(int));
        assert(f);
        f[0] = eval[0];
        return f;
    }
    int *eval_r0;
    int *eval_r1;
    array_split(&eval_r0, &eval_r1, eval, d);
    int *r0 = inv_fft(eval_r0, d / 2, zp_mul(omega, omega));
    int *r1 = inv_fft(eval_r1, d / 2, zp_mul(omega, omega));
    free(eval_r0);
    free(eval_r1);
    int omega_i = omega;
    for (int i = 1; i <= d / 2; i++)
    {
        r1[d / 2 - i] = zp_opp(zp_mul(r1[d / 2 - i], omega_i));
        omega_i = zp_mul(omega_i, omega);
    }
    int *f = (int *)malloc(d * sizeof(int));
    assert(f);
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
    int q, d;
    for (q = p - 1, d = 1; (q & 1) == 0; q >>= 1, d <<= 1)
        ;
    int omega = zp_prim_root_min(q, d);
    int *f_coeffs = inv_fft(eval, d, omega);
    int degree = d - 1;
    while (degree >= 0 && f_coeffs[degree] == 0)
        degree--;
    int *coeffs;
    if (degree == -1)
    {
        coeffs = (int *)calloc(1, sizeof(int));
        assert(coeffs);
    }
    else
    {
        coeffs = (int *)malloc((degree + 1) * sizeof(int));
        assert(coeffs);
        for (int i = 0; i <= degree; i++)
            coeffs[i] = f_coeffs[i];
    }
    free(f_coeffs);
    poly_set(f, degree, coeffs);
}

void poly_mul_fft(poly *rop, poly *op1, poly *op2)
{
    int q, d;
    for (q = p - 1, d = 1; (q & 1) == 0; q >>= 1, d <<= 1)
        ;
    int *eval_op1 = NULL;
    poly_fft(op1, &eval_op1);
    int *eval_op2 = NULL;
    poly_fft(op2, &eval_op2);
    int *eval_rop = (int *)malloc(d * sizeof(int));
    assert(eval_rop);
    for (int i = 0; i < d; i++)
        eval_rop[i] = zp_mul(eval_op1[i], eval_op2[i]);
    poly_inv_fft(rop, eval_rop);
    free(eval_op1);
    free(eval_op2);
    free(eval_rop);
}