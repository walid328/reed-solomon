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

void set_poly(poly *f, int deg, int *coeffs)
{
    f->degree = deg;
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

char *str_poly_min(poly *q)
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
        if (q->coeffs[i] != 0)
        {
            if (index != 0)
            {
                poly_string[index + 0] = ' ';
                poly_string[index + 1] = '+';
                poly_string[index + 2] = ' ';
                index += 3;
            }
            str_add_int(poly_string, &index, q->coeffs[i]);
            if (i > 0)
            {
                poly_string[index + 0] = '*';
                poly_string[index + 1] = 'x';
                index += 2;
                if (i > 1)
                {
                    poly_string[index + 0] = '^';
                    index += 1;
                    str_add_int(poly_string, &index, i);
                }
            }
        }
    }
    poly_string[index] = '\0';
    return poly_string;
}

void print_poly_iso_length(poly *q)
{
    char *poly_string = str_poly_iso_length(q);
    printf("%s\n", poly_string);
    free(poly_string);
}

void print_poly_min(poly *q)
{
    char *poly_string = str_poly_min(q);
    printf("%s\n", poly_string);
    free(poly_string);
}

void print_poly(poly *q)
{
    if (q->degree == -1)
        printf("0");
    for (int i = 0; i <= q->degree; i++)
    {
        if (i != 0)
            printf(" + ");
        printf("%d", q->coeffs[i]);
        if (i == 1)
            printf("*x");
        if (i > 1)
            printf("*x^%d", i);
    }
    printf("\n");
}

/******************************************************/

poly *new_poly(void)
{
    poly *f = (poly *)malloc(sizeof(poly));
    assert(f);
    set_poly(f, -1, NULL);
    return f;
}

poly *new_poly_0(void)
{
    poly *f = new_poly();
    int *coeffs = (int *)calloc(1, sizeof(int));
    assert(coeffs);
    set_poly(f, -1, coeffs);
    return f;
}

poly *new_poly_1(void)
{
    poly *f = new_poly();
    int *coeffs = (int *)malloc(1 * sizeof(int));
    assert(coeffs);
    coeffs[0] = 1;
    set_poly(f, 0, coeffs);
    return f;
}

poly *new_poly_from_coeffs(int deg, ...)
{
    if (deg == -1)
        return new_poly_0();
    poly *f = new_poly();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    assert(coeffs);
    va_list valist;
    va_start(valist, deg);
    for (int i = 0; i <= deg; i++)
        coeffs[i] = va_arg(valist, int);
    va_end(valist);
    set_poly(f, deg, coeffs);
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

poly *new_poly_from_str(char *str)
{
    int deg = -1;
    int *coeffs_str = coeffs_from_str(str, &deg);
    if (deg == -1)
    {
        free(coeffs_str);
        return new_poly_0();
    }
    poly *f = new_poly();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    assert(coeffs);
    for (int i = 0; i <= deg; i++)
        coeffs[i] = coeffs_str[i];
    free(coeffs_str);
    set_poly(f, deg, coeffs);
    return f;
}

poly *new_poly_from_copy(poly *source)
{
    int deg = source->degree;
    if (deg == -1)
        return new_poly_0();
    poly *copy = new_poly();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    assert(coeffs);
    for (int i = 0; i <= deg; i++)
        coeffs[i] = source->coeffs[i];
    set_poly(copy, deg, coeffs);
    return copy;
}

poly *new_rand_poly(int deg)
{
    poly *f = new_poly();
    int *coeffs = (int *)malloc((deg + 1) * sizeof(int));
    assert(coeffs);
    for (int i = 0; i <= deg; i++)
        coeffs[i] = rand_zp();
    while (coeffs[deg] == 0)
        coeffs[deg] = rand_zp();
    set_poly(f, deg, coeffs);
    return f;
}

/******************************************************/

void clear_poly(poly *f)
{
    f->degree = -1;
    f->coeffs = NULL;
}

void clear_full_poly(poly *f)
{
    f->degree = -1;
    assert(f->coeffs);
    free(f->coeffs);
    f->coeffs = NULL;
}

void free_poly(poly *f)
{
    assert(f);
    free(f);
}

void free_full_poly(poly *f)
{
    assert(f);
    assert(f->coeffs);
    free(f->coeffs);
    free(f);
}

/******************************************************/

void add_poly(poly *rop, poly *op1, poly *op2)
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
            coeffs[i] = add_zp(op1->coeffs[i], op2->coeffs[i]);
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
        for (i = op1->degree; i >= 0 && add_zp(op1->coeffs[i], op2->coeffs[i]) == 0; i--)
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
                coeffs[j] = add_zp(op1->coeffs[j], op2->coeffs[j]);
        }
    }
    if (rop == op1)
        clear_full_poly(op1);
    if (rop == op2)
        clear_full_poly(op2);
    set_poly(rop, deg, coeffs);
}

void sub_poly(poly *rop, poly *op1, poly *op2)
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
            coeffs[i] = sub_zp(op1->coeffs[i], op2->coeffs[i]);
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
            coeffs[i] = sub_zp(op1->coeffs[i], op2->coeffs[i]);
        for (; i <= deg; i++)
            coeffs[i] = opp_zp(op2->coeffs[i]);
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
        for (i = op1->degree; i >= 0 && sub_zp(op1->coeffs[i], op2->coeffs[i]) == 0; i--)
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
                coeffs[j] = sub_zp(op1->coeffs[j], op2->coeffs[j]);
        }
    }
    if (rop == op1)
        clear_full_poly(op1);
    if (rop == op2)
        clear_full_poly(op2);
    set_poly(rop, deg, coeffs);
}

void mul_poly(poly *rop, poly *op1, poly *op2)
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
                coeffs[i + j] = add_zp(coeffs[i + j], mul_zp(op1->coeffs[i], op2->coeffs[j]));
    }
    if (rop == op1)
        clear_full_poly(op1);
    if (rop == op2)
        clear_full_poly(op2);
    set_poly(rop, deg, coeffs);
}

void mul_scalar_poly(poly *rop, int op1, poly *op2)
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
            coeffs[i] = mul_zp(op2->coeffs[i], op1);
    }
    if (rop == op2)
        clear_full_poly(op2);
    set_poly(rop, deg, coeffs);
}

void euclid_div_poly(poly *q, poly *r, poly *op1, poly *op2)
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
    int inv_LC_op2 = inv_zp(op2->coeffs[op2->degree]);
    while (deg_r >= op2->degree)
    {
        int c = mul_zp(rest[deg_r], inv_LC_op2);
        coeffs_q[deg_r - op2->degree] = c;
        for (int i = 0; i <= op2->degree; i++)
            rest[i + deg_r - op2->degree] = sub_zp(rest[i + deg_r - op2->degree], mul_zp(c, op2->coeffs[i]));
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
    set_poly(q, deg_q, coeffs_q);
    set_poly(r, deg_r, coeffs_r);
}

void xgcd_poly(poly **d, poly **u, poly **v, poly *op1, poly *op2)
{
    *d = new_poly_from_copy(op1);
    *u = new_poly_1();
    *v = new_poly_0();
    poly *r1 = new_poly_from_copy(op2);
    poly *u1 = new_poly_0();
    poly *v1 = new_poly_1();
    while (r1->degree > -1)
    {
        poly *q = new_poly();
        poly *r2 = new_poly();
        euclid_div_poly(q, r2, *d, r1);
        poly *u2 = new_poly();
        mul_poly(u2, q, u1);
        sub_poly(u2, *u, u2);
        poly *v2 = new_poly();
        mul_poly(v2, q, v1);
        sub_poly(v2, *v, v2);
        free_full_poly(q);
        free_full_poly(*d);
        free_full_poly(*u);
        free_full_poly(*v);
        *d = r1;
        *u = u1;
        *v = v1;
        r1 = r2;
        u1 = u2;
        v1 = v2;
    }
    int inv_LC_d = inv_zp((*d)->coeffs[(*d)->degree]);
    mul_scalar_poly(*d, inv_LC_d, *d);
    mul_scalar_poly(*u, inv_LC_d, *u);
    mul_scalar_poly(*v, inv_LC_d, *v);
    free_full_poly(r1);
    free_full_poly(u1);
    free_full_poly(v1);
}

/******************************************************/

void deriv_poly(poly *rop, poly *op)
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
            coeffs[i] = mul_zp(i + 1, op->coeffs[i + 1]);
    }
    if (rop == op)
        clear_full_poly(op);
    set_poly(rop, deg, coeffs);
}

int eval_poly(poly *f, int x)
{
    int y = 0;
    int x_i = 1;
    for (int i = 0; i <= f->degree; i++)
    {
        y = add_zp(y, mul_zp(x_i, f->coeffs[i]));
        x_i = mul_zp(x_i, x);
    }
    return y;
}

void multi_eval_poly(poly *f, int n, int *a, int *b)
{
    for (int i = 0; i <= n; i++)
        b[i] = eval_poly(f, a[i]);
}

poly *interpolation(int *a, int *b, int n)
{
    poly *f = new_poly_1();
    for (int i = 0; i <= n; i++)
    {
        poly *x_m_ai = new_poly_from_coeffs(1, opp_zp(a[i]), 1);
        mul_poly(f, f, x_m_ai);
        free_full_poly(x_m_ai);
    }
    poly *fd = new_poly();
    deriv_poly(fd, f);
    poly *g = new_poly_0();
    for (int i = 0; i <= n; i++)
    {
        poly *x_m_ai = new_poly_from_coeffs(1, opp_zp(a[i]), 1);
        int fd_ai = eval_poly(fd, a[i]);
        mul_scalar_poly(x_m_ai, fd_ai, x_m_ai);
        poly *l_i = new_poly();
        poly *r = new_poly();
        euclid_div_poly(l_i, r, f, x_m_ai);
        mul_scalar_poly(l_i, b[i], l_i);
        add_poly(g, g, l_i);
        free_full_poly(x_m_ai);
        free_full_poly(l_i);
        free_full_poly(r);
    }
    free_full_poly(f);
    free_full_poly(fd);
    return g;
}