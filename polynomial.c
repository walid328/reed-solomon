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

int poly_leading_coeff(const poly f)
{
    if (f->deg == -1)
        return 0;
    return f->coeffs[f->deg];
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
    if (poly_is_zero(f))
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
    poly f = (poly)malloc(sizeof(struct polynomial));
    assert(f);
    poly_set(f, -1, NULL);
    return f;
}

poly poly_new_deg(int deg)
{
    poly f = poly_new();
    poly_set_deg(f, deg);
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

poly poly_new_rand(int deg)
{
    poly f = poly_new();
    if (deg >= 0)
    {
        array coeffs = array_new(deg + 1);
        for (int i = 0; i <= deg; i++)
            coeffs[i] = zp_rand();
        while (coeffs[deg] == 0)
            coeffs[deg] = zp_rand();
        poly_set(f, deg, coeffs);
    }
	return f;
}

poly poly_new_copy(const poly src)
{
    poly f = poly_new();
    poly_copy(f, src);
    return f;
}

/******************************************************/

void poly_clear(poly f)
{
    array_free(f->coeffs);
    poly_set(f, -1, NULL);
}

void poly_clear_multi(int qty, ...)
{
    va_list valist;
    va_start(valist, qty);
    for (int i = 0; i < qty; i++)
    {
        poly f = va_arg(valist, poly);
        poly_clear(f);
    }
    va_end(valist);
}

void poly_free(poly f)
{
    poly_clear(f);
    free(f);
}

void poly_free_multi(int qty, ...)
{
    va_list valist;
    va_start(valist, qty);
    for (int i = 0; i < qty; i++)
    {
        poly f = va_arg(valist, poly);
        poly_free(f);
    }
    va_end(valist);
}

/******************************************************/

bool poly_is_zero(const poly f)
{
    return f->deg == -1;
}

bool poly_equal(const poly f, const poly g)
{
    if (poly_deg(f) != poly_deg(g))
        return false;
    for (int i = 0; i < poly_deg(f); i++)
        if (poly_coeffs(f)[i] != poly_coeffs(g)[i])
            return false;
    return true;
}

void poly_copy(poly dst, const poly src)
{
    int deg = src->deg;
    array coeffs = NULL;
    if (deg > -1)
    {
        coeffs = array_new(deg + 1);
        memcpy(coeffs, src->coeffs, (deg + 1) * sizeof(int));
    }
    poly_clear(dst);
    poly_set(dst, deg, coeffs);
}

void poly_set_deg(poly f, int deg)
{
    poly_clear(f);
    if (deg > -1)
    {
        f->deg = deg;
        f->coeffs = array_new(deg + 1);
    }
}

void poly_set_coeffs(poly f, int deg, ...)
{
    poly_set_deg(f, deg);
    if (deg > -1)
    {
        va_list valist;
        va_start(valist, deg);
        for (int i = 0; i <= deg; i++)
            f->coeffs[i] = va_arg(valist, int);
        va_end(valist);
    }
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
    else if (poly_is_zero(op1))
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
    else if (poly_is_zero(op1))
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
    if (poly_is_zero(op1) || poly_is_zero(op2))
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
    if (op1 == 0 || poly_is_zero(op2))
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
    poly_clear(rop);
    poly_set(rop, deg, coeffs);
}

void poly_euc_div(poly q, poly r, const poly op1, const poly op2)
{
    if (op1->deg < op2->deg)
    {
        poly_clear(q);
        poly_copy(r, op1);
        return;
    }
    if (poly_is_zero(op2))
    {
        fprintf(stderr, "Can't divide by 0!\n");
        exit(EXIT_FAILURE);
    }
    int deg_q, *coeffs_q;
    int deg_r, *coeffs_r;
    array rest = array_new(op1->deg + 1);
    for (int i = 0; i <= op1->deg; i++)
        rest[i] = op1->coeffs[i];
    deg_q = op1->deg - op2->deg;
    coeffs_q = array_new_zeros(deg_q + 1);
    deg_r = op1->deg;
    int inv_LC_op2 = zp_inv(poly_leading_coeff(op2));
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
    poly_clear(q);
    poly_clear(r);
    poly_set(q, deg_q, coeffs_q);
    poly_set(r, deg_r, coeffs_r);
}

void poly_xgcd(poly d, poly u, poly v, const poly op1, const poly op2)
{
    poly r0 = poly_new_copy(op1);
    poly u0 = poly_new_set(0, 1);
    poly v0 = poly_new();
    poly r1 = poly_new_copy(op2);
    poly u1 = poly_new();
    poly v1 = poly_new_set(0, 1);
    poly q = poly_new();
    while (r1->deg > -1)
    {
        poly r2 = poly_new();
        poly u2 = poly_new();
        poly v2 = poly_new();
        poly_euc_div(q, r2, r0, r1);
        poly_mul(u2, q, u1);
        poly_sub(u2, u0, u2);
        poly_mul(v2, q, v1);
        poly_sub(v2, v0, v2);
        poly_free_multi(3, r0, u0, v0);
        r0 = r1;
        u0 = u1;
        v0 = v1;
        r1 = r2;
        u1 = u2;
        v1 = v2;
    }
    poly_copy(d, r0);
    poly_copy(u, u0);
    poly_copy(v, v0);
    poly_free_multi(7, r0, u0, v0, r1, u1, v1, q);
}

void poly_xgcd_partial(poly d, poly u, poly v, const poly op1, const poly op2, int limit)
{
    poly r0 = poly_new_copy(op1);
    poly u0 = poly_new_set(0, 1);
    poly v0 = poly_new();
    poly r1 = poly_new_copy(op2);
    poly u1 = poly_new();
    poly v1 = poly_new_set(0, 1);
    poly q = poly_new();
    while (r1->deg >= limit)
    {
        poly r2 = poly_new();
        poly u2 = poly_new();
        poly v2 = poly_new();
        poly_euc_div(q, r2, r0, r1);
        poly_mul(u2, q, u1);
        poly_sub(u2, u0, u2);
        poly_mul(v2, q, v1);
        poly_sub(v2, v0, v2);
        poly_free_multi(3, r0, u0, v0);
        r0 = r1;
        u0 = u1;
        v0 = v1;
        r1 = r2;
        u1 = u2;
        v1 = v2;
    }
    poly_copy(d, r1);
    poly_copy(u, u1);
    poly_copy(v, v1);
    poly_free_multi(7, r0, u0, v0, r1, u1, v1, q);
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
        poly_free_multi(3, x_m_ai, l_i, r);
    }
    poly_free_multi(2, f, fd);
}

array poly_dft(const poly f, int d)
{
	if (d > n && n % d != 0)
	{
		fprintf(stderr, "There is no %d-th root of unity!\n", d);
		exit(EXIT_FAILURE);
	}
    array eval = array_new(d);
    for (int i = 0; i < d; i++)
        eval[i] = poly_eval(f, omegas[i * (n / d)]);
    return eval;
}

void poly_inv_dft(poly rop, array eval, int d)
{
	if (d > n && n % d != 0)
	{
		fprintf(stderr, "There is no %d-th root of unity!\n", d);
		exit(EXIT_FAILURE);
	}
    array points = array_new(d);
    for (int i = 0; i < d; i++)
        points[i] = omegas[i * (n / d)];
    interpolation(rop, points, eval, d);
    array_free(points);
}

/******************************************************/

array poly_fft(const poly f, int d)
{
	if (d > n && n % d != 0)
	{
		fprintf(stderr, "There is no %d-th root of unity!\n", d);
		exit(EXIT_FAILURE);
	}
    // We need some padding.
    array coeffs = array_new_zeros(d);
    for (int i = 0; i <= f->deg; i++)
        coeffs[i] = f->coeffs[i];
    array eval = array_fft(coeffs, d);
    array_free(coeffs);
    return eval;
}

void poly_inv_fft(poly rop, array eval, int d)
{
	if (d > n && n % d != 0)
	{
		fprintf(stderr, "There is no %d-th root of unity!\n", d);
		exit(EXIT_FAILURE);
	}
    array f_coeffs = array_inv_fft(eval, d);
    int deg = d - 1;
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
    poly_clear(rop);
    poly_set(rop, deg, coeffs);
}

/******************************************************/

bool fast_mul_usefull(int deg1, int deg2)
{
	if (deg1 < 2 || deg2 < 2)
		return false;
	long long int dn = deg1;
    long long int dm = deg2;
    long long int logn = 1;
    long long int dlogn = 2;
    while (dlogn < dn)
    {
        logn++;
        dlogn <<= 1;
    }
	return (dn * logn) < (dn * dm);
}

void poly_fast_mul(poly rop, const poly op1, const poly op2)
{
    if (poly_is_zero(op1) || poly_is_zero(op2))
    {
        poly_clear(rop);
    }
	else if (!fast_mul_usefull(op1->deg, op2->deg))
	{
		poly_mul(rop, op1, op2);
	}
    else
    {
        int d = 1;
        while (d < op1->deg + op2->deg + 1)
            d *= 2;
        if (d == 2 * n)
        {
			int d11 = op1->deg / 2;
			int d12 = (op1->deg - 1) / 2;
			int d21 = op2->deg / 2;
			int d22 = (op2->deg - 1) / 2;
			int x11 = d11 + 1;
			int x21 = d21 + 1;
			while (d11 >= 0 && op1->coeffs[d11] == 0)
				d11--;
			while (d21 >= 0 && op2->coeffs[d21] == 0)
				d21--;
            poly op11 = poly_new_deg(d11);
			poly op12 = poly_new_deg(d12);
            poly op21 = poly_new_deg(d21);
			poly op22 = poly_new_deg(d22);
			for (int i = 0; i <= d11; i++)
				op11->coeffs[i] = op1->coeffs[i];
			for (int i = x11; i - x11 <= d12; i++)
				op12->coeffs[i - x11] = op1->coeffs[i];
			for (int i = 0; i <= d21; i++)
				op21->coeffs[i] = op2->coeffs[i];
			for (int i = x21; i - x21 <= d22; i++)
				op22->coeffs[i - x21] = op2->coeffs[i];
			poly op11op21 = poly_new();
			poly op11op22 = poly_new();
			poly op12op21 = poly_new();
			poly op12op22 = poly_new();
			poly_fast_mul(op11op21, op11, op21);
			poly_fast_mul(op11op22, op11, op22);
			poly_fast_mul(op12op21, op12, op21);
			poly_fast_mul(op12op22, op12, op22);
			poly_set_deg(rop, op1->deg + op2->deg);
			for (int i = 0; i <= op1->deg + op2->deg; i++)
			{
				int coeff = 0;
				if (i <= d11 + d21)
					coeff = zp_add(coeff, op11op21->coeffs[i]);
				if (x21 <= i && i - x21 <= d11 + d22)
					coeff = zp_add(coeff, op11op22->coeffs[i - x21]);
				if (x11 <= i && i - x11 <= d12 + d21)
					coeff = zp_add(coeff, op12op21->coeffs[i - x11]);
				if (x11 + x21 <= i && i - x11 - x21 <= d12 + d22)
					coeff = zp_add(coeff, op12op22->coeffs[i - x11 - x21]);
				rop->coeffs[i] = coeff;
			}
			poly_free_multi(8, op11, op12, op21, op22, op11op21, op11op22, op12op21, op12op22);
		}
		else
		{
			array eval_op1 = poly_fft(op1, d);
			array eval_op2 = poly_fft(op2, d);
			array eval_rop = array_new(d);
			for (int i = 0; i < d; i++)
				eval_rop[i] = zp_mul(eval_op1[i], eval_op2[i]);
			poly_inv_fft(rop, eval_rop, d);
			array_free(eval_op1);
			array_free(eval_op2);
			array_free(eval_rop);
		}
    }
}

bool fast_euc_div_useful(int deg1, int deg2)
{
    long long int dn = deg1;
    long long int dm = deg2;
    long long int logn = 1;
    long long int dlogn = 2;
    while (dlogn < dn)
    {
        logn++;
        dlogn <<= 1;
    }
    long long int loglogn = 1;
    long long int dloglogn = 2;
    while (dloglogn < logn)
    {
        loglogn++;
        dloglogn <<= 1;
    }
    return (dn * logn * loglogn) < (dm * (dn - dm));
}

void poly_fast_euc_div(poly quo, poly rem, const poly op1, const poly op2)
{
    if (op1->deg < op2->deg)
    {
        poly_clear(quo);
        poly_copy(rem, op1);
        return;
    }
    if (poly_is_zero(op2))
    {
        fprintf(stderr, "Can't divide by 0!\n");
        exit(EXIT_FAILURE);
    }
    if (fast_euc_div_useful(op1->deg, op2->deg))
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
	else
        poly_euc_div(quo, rem, op1, op2);
}

poly *poly_half_gcd(const poly r0, const poly r1)
{
    // r0 = R0, r1 = R1
    // 1.
    int m = (r0->deg + 1) / 2;
    // 2.
    if (r1->deg < m)
    {
        poly *mdpgcd = (poly *)malloc(4 * sizeof(poly));
        assert(mdpgcd);
        mdpgcd[0] = poly_new_set(0, 1);
        mdpgcd[1] = poly_new();
        mdpgcd[2] = poly_new();
        mdpgcd[3] = poly_new_set(0, 1);
        // 3.
        return mdpgcd;
    }
    // 4.
    // r0s = R0*
    poly r0s = poly_new_deg(r0->deg - m);
    for (int i = m; i <= r0->deg; i++)
        r0s->coeffs[i - m] = r0->coeffs[i];
    // r1s = R1*
    poly r1s = poly_new_deg(r1->deg - m);
    for (int i = m; i <= r1->deg; i++)
        r1s->coeffs[i - m] = r1->coeffs[i];
    // 5.
    // d0s = D0*
    poly *d0s = poly_half_gcd(r0s, r1s);
    // 6.
    poly prod1 = poly_new();
    poly prod2 = poly_new();
    poly_fast_mul(prod1, d0s[0], r0);
    poly_fast_mul(prod2, d0s[1], r1);
    // rj = Rj
    poly rj = poly_new();
    poly_add(rj, prod1, prod2);
    poly_fast_mul(prod1, d0s[2], r0);
    poly_fast_mul(prod2, d0s[3], r1);
    // rjp1 = R(j+1)
    poly rjp1 = poly_new();
    poly_add(rjp1, prod1, prod2);
    // 7.
    if (rjp1->deg < m)
    {
        poly_free_multi(6, r0s, r1s, prod1, prod2, rj, rjp1);
        // 8.
        return d0s;
    }
    // 9.
    // qj = Qj
    poly qj = poly_new();
    // rjp2 = R(j+2)
    poly rjp2 = poly_new();
    poly_fast_euc_div(qj, rjp2, rj, rjp1);
    // Tj = [[0, 1], [1, -qj]]
    poly w2 = poly_new();
    poly_fast_mul(prod1, qj, d0s[2]);
    poly_sub(w2, d0s[0], prod1);
    poly w3 = poly_new();
    poly_fast_mul(prod1, qj, d0s[3]);
    poly_sub(w3, d0s[1], prod1);
    // Tj * D0* = [[d0s[2], d0s[3]], [w2, w3]]
    poly *mdpgcd = (poly *)malloc(4 * sizeof(poly));
    assert(mdpgcd);
    // 10.
    if (rjp2->deg < m)
    {
        mdpgcd[0] = d0s[2];
        mdpgcd[1] = d0s[3];
        mdpgcd[2] = w2;
        mdpgcd[3] = w3;
        poly_free_multi(10, r0s, r1s, prod1, prod2, rj, rjp1, qj, rjp2, d0s[0], d0s[1]);
        free(d0s);
        // 11.
        return mdpgcd;
    }
    // 12.
    int l = 2 * m - rjp1->deg;
    // 13.
    // rjp1s = R(j+1)*
    poly rjp1s = poly_new_deg(rjp1->deg - l);
    for (int i = l; i <= rjp1->deg; i++)
        rjp1s->coeffs[i - l] = rjp1->coeffs[i];
    // rjp2s = R(j+2)*
    poly rjp2s;
    if (rjp2->deg < l)
        rjp2s = poly_new();
    else
    {
        rjp2s = poly_new_deg(rjp2->deg - l);
        for (int i = l; i <= rjp2->deg; i++)
            rjp2s->coeffs[i - l] = rjp2->coeffs[i];
    }
    // 14.
    // d1s = D1*
    poly *d1s = poly_half_gcd(rjp1s, rjp2s);
    // 15.
    for (int i = 0; i < 4; i++)
        mdpgcd[i] = poly_new();
    poly_fast_mul(prod1, d1s[0], d0s[2]);
    poly_fast_mul(prod2, d1s[1], w2);
    poly_add(mdpgcd[0], prod1, prod2);
    poly_fast_mul(prod1, d1s[0], d0s[3]);
    poly_fast_mul(prod2, d1s[1], w3);
    poly_add(mdpgcd[1], prod1, prod2);
    poly_fast_mul(prod1, d1s[2], d0s[2]);
    poly_fast_mul(prod2, d1s[3], w2);
    poly_add(mdpgcd[2], prod1, prod2);
    poly_fast_mul(prod1, d1s[2], d0s[3]);
    poly_fast_mul(prod2, d1s[3], w3);
    poly_add(mdpgcd[3], prod1, prod2);
    for (int i = 0; i < 4; i++)
        poly_free_multi(2, d0s[i], d1s[i]);
    free(d0s);
    free(d1s);
    poly_free_multi(12, r0s, r1s, prod1, prod2, rj, rjp1, qj, rjp2, w2, w3, rjp1s, rjp2s);
    return mdpgcd;
}

poly *poly_fast_gcd_matrix(const poly r0, const poly r1)
{
    // r0 = R0, r1 = R1
    // 1.
    // mdpgcd = D
    poly *mdpgcd = poly_half_gcd(r0, r1);
    // 2.
    poly prod1 = poly_new();
    poly prod2 = poly_new();
    poly_fast_mul(prod1, mdpgcd[0], r0);
    poly_fast_mul(prod2, mdpgcd[1], r1);
    // rj = Rj
    poly rj = poly_new();
    poly_add(rj, prod1, prod2);
    poly_fast_mul(prod1, mdpgcd[2], r0);
    poly_fast_mul(prod2, mdpgcd[3], r1);
    // rjp1 = R(j+1)
    poly rjp1 = poly_new();
    poly_add(rjp1, prod1, prod2);
    // 3.
    if (poly_is_zero(rjp1))
    {
        poly_free_multi(4, prod1, prod2, rj, rjp1);
        return mdpgcd;
    }
    // 4.
    // qj = Qj
    poly qj = poly_new();
    // rjp2 = R(j+2)
    poly rjp2 = poly_new();
    poly_fast_euc_div(qj, rjp2, rj, rjp1);
    // 5.
    // Tj = [[0, 1], [1, -qj]]
    // 6.
    poly w2 = poly_new();
    poly_fast_mul(prod1, qj, mdpgcd[2]);
    poly_sub(w2, mdpgcd[0], prod1);
    poly w3 = poly_new();
    poly_fast_mul(prod1, qj, mdpgcd[3]);
    poly_sub(w3, mdpgcd[1], prod1);
    // Tj * D = [[mdpgcd[2], mdpgcd[3]], [w2, w3]]
    // mgcd = matrix of gcd of R0 and R1
    poly *mgcd = (poly *)malloc(4 * sizeof(poly));
    assert(mgcd);
    if (poly_is_zero(rjp2))
    {
        mgcd[0] = mdpgcd[2];
        mgcd[1] = mdpgcd[3];
        mgcd[2] = w2;
        mgcd[3] = w3;
        poly_free_multi(8, mdpgcd[0], mdpgcd[1], prod1, prod2, rj, rjp1, rjp2, qj);
        free(mdpgcd);
        return mgcd;
    }
    // 7.
    // mrjp12 = N
    poly *mrjp12 = poly_fast_gcd_matrix(rjp1, rjp2);
    // 8.
    for (int i = 0; i < 4; i++)
        mgcd[i] = poly_new();
    poly_fast_mul(prod1, mrjp12[0], mdpgcd[2]);
    poly_fast_mul(prod2, mrjp12[1], w2);
    poly_add(mgcd[0], prod1, prod2);
    poly_fast_mul(prod1, mrjp12[0], mdpgcd[3]);
    poly_fast_mul(prod2, mrjp12[1], w3);
    poly_add(mgcd[1], prod1, prod2);
    poly_fast_mul(prod1, mrjp12[2], mdpgcd[2]);
    poly_fast_mul(prod2, mrjp12[3], w2);
    poly_add(mgcd[2], prod1, prod2);
    poly_fast_mul(prod1, mrjp12[2], mdpgcd[3]);
    poly_fast_mul(prod2, mrjp12[3], w3);
    poly_add(mgcd[3], prod1, prod2);
    for (int i = 0; i < 4; i++)
        poly_free_multi(2, mdpgcd[i], mrjp12[i]);
    free(mdpgcd);
    free(mrjp12);
    poly_free_multi(8, prod1, prod2, rj, rjp1, rjp2, qj, w2, w3);
    return mgcd;
}

void poly_fast_xgcd(poly d, poly u, poly v, const poly op1, const poly op2)
{
    if (poly_is_zero(op1) && poly_is_zero(op2))
    {
        fprintf(stderr, "there is no gcd of 0 and 0\n");
        exit(EXIT_FAILURE);
    }
    if (op1->deg == op2->deg)
    {
        poly quo = poly_new();
        poly rem = poly_new();
        poly_fast_euc_div(quo, rem, op2, op1);
        poly_fast_xgcd(d, u, v, op1, rem);
        poly qv = poly_new();
        poly_fast_mul(qv, quo, v);
        poly_sub(u, u, qv);
        poly_free_multi(3, quo, rem, qv);
    }
    else
    {
        if (op1->deg > op2->deg)
        {
            poly *mab = poly_fast_gcd_matrix(op1, op2);
            poly prod1 = poly_new();
            poly prod2 = poly_new();
            poly_fast_mul(prod1, mab[0], op1);
            poly_fast_mul(prod2, mab[1], op2);
            poly_add(d, prod1, prod2);
            poly_copy(u, mab[0]);
            poly_copy(v, mab[1]);
            poly_free_multi(6, mab[0], mab[1], mab[2], mab[3], prod1, prod2);
            free(mab);
        }
        else if (op1->deg < op2->deg)
        {
            poly *mab = poly_fast_gcd_matrix(op2, op1);
            poly prod1 = poly_new();
            poly prod2 = poly_new();
            poly_fast_mul(prod1, mab[0], op2);
            poly_fast_mul(prod2, mab[1], op1);
            poly_add(d, prod1, prod2);
            poly_copy(u, mab[1]);
            poly_copy(v, mab[0]);
            poly_free_multi(6, mab[0], mab[1], mab[2], mab[3], prod1, prod2);
            free(mab);
        }
    }
}

poly *poly_fast_gcd_partial_matrix(const poly r0, const poly r1, int limit)
{
    if ((r0->deg + 1) / 2 > limit)
    {
        // r0 = R0, r1 = R1
        // 1.
        // mdpgcd = D
        poly *mdpgcd = poly_half_gcd(r0, r1);
        // 2.
        poly prod1 = poly_new();
        poly prod2 = poly_new();
        poly_fast_mul(prod1, mdpgcd[0], r0);
        poly_fast_mul(prod2, mdpgcd[1], r1);
        // rj = Rj
        poly rj = poly_new();
        poly_add(rj, prod1, prod2);
        poly_fast_mul(prod1, mdpgcd[2], r0);
        poly_fast_mul(prod2, mdpgcd[3], r1);
        // rjp1 = R(j+1)
        poly rjp1 = poly_new();
        poly_add(rjp1, prod1, prod2);
        // 3.
        if (rjp1->deg < limit)
        {
            poly_free_multi(4, prod1, prod2, rj, rjp1);
            return mdpgcd;
        }
        // 4.
        // qj = Qj
        poly qj = poly_new();
        // rjp2 = R(j+2)
        poly rjp2 = poly_new();
        poly_fast_euc_div(qj, rjp2, rj, rjp1);
        // 5.
        // Tj = [[0, 1], [1, -qj]]
        // 6.
        poly w2 = poly_new();
        poly_fast_mul(prod1, qj, mdpgcd[2]);
        poly_sub(w2, mdpgcd[0], prod1);
        poly w3 = poly_new();
        poly_fast_mul(prod1, qj, mdpgcd[3]);
        poly_sub(w3, mdpgcd[1], prod1);
        // Tj * D = [[mdpgcd[2], mdpgcd[3]], [w2, w3]]
        // mgcd = matrix of gcd of R0 and R1
        poly *mgcd = (poly *)malloc(4 * sizeof(poly));
        assert(mgcd);
        if (rjp2->deg < limit)
        {
            mgcd[0] = mdpgcd[2];
            mgcd[1] = mdpgcd[3];
            mgcd[2] = w2;
            mgcd[3] = w3;
            poly_free_multi(8, mdpgcd[0], mdpgcd[1], prod1, prod2, rj, rjp1, rjp2, qj);
            free(mdpgcd);
            return mgcd;
        }
        // 7.
        // mrjp12 = N
        poly *mrjp12 = poly_fast_gcd_partial_matrix(rjp1, rjp2, limit);
        // 8.
        for (int i = 0; i < 4; i++)
            mgcd[i] = poly_new();
        poly_fast_mul(prod1, mrjp12[0], mdpgcd[2]);
        poly_fast_mul(prod2, mrjp12[1], w2);
        poly_add(mgcd[0], prod1, prod2);
        poly_fast_mul(prod1, mrjp12[0], mdpgcd[3]);
        poly_fast_mul(prod2, mrjp12[1], w3);
        poly_add(mgcd[1], prod1, prod2);
        poly_fast_mul(prod1, mrjp12[2], mdpgcd[2]);
        poly_fast_mul(prod2, mrjp12[3], w2);
        poly_add(mgcd[2], prod1, prod2);
        poly_fast_mul(prod1, mrjp12[2], mdpgcd[3]);
        poly_fast_mul(prod2, mrjp12[3], w3);
        poly_add(mgcd[3], prod1, prod2);
        for (int i = 0; i < 4; i++)
            poly_free_multi(2, mdpgcd[i], mrjp12[i]);
        free(mdpgcd);
        free(mrjp12);
        poly_free_multi(8, prod1, prod2, rj, rjp1, rjp2, qj, w2, w3);
        return mgcd;
    }
    else
    {
        int t = 2 * limit - r0->deg;
        poly r0s = poly_new_deg(r0->deg - t);
        for (int i = t; i <= r0->deg; i++)
            r0s->coeffs[i - t] = r0->coeffs[i];
        poly r1s = poly_new_deg(r1->deg - t);
        for (int i = t; i <= r1->deg; i++)
            r1s->coeffs[i - t] = r1->coeffs[i];
        poly *ds = poly_half_gcd(r0s, r1s);
        poly_free_multi(2, r0s, r1s);
        return ds;
    }
}

void poly_fast_xgcd_partial(poly d, poly u, poly v, const poly op1, const poly op2, int limit)
{
    if (poly_is_zero(op1) && poly_is_zero(op2))
    {
        fprintf(stderr, "there is no gcd of 0 and 0\n");
        exit(EXIT_FAILURE);
    }
    if (op1->deg < limit || op2->deg < limit)
    {
        if (op1->deg >= op2->deg)
        {
            poly_copy(d, op1);
            poly_set_coeffs(u, 0, 1);
            poly_clear(v);
        }
        else
        {
            poly_copy(d, op2);
            poly_clear(u);
            poly_set_coeffs(v, 0, 1);
        }
    }
    if (op1->deg == op2->deg)
    {
        poly quo = poly_new();
        poly rem = poly_new();
        poly_fast_euc_div(quo, rem, op2, op1);
        poly_fast_xgcd_partial(d, u, v, op1, rem, limit);
        poly qv = poly_new();
        poly_fast_mul(qv, quo, v);
        poly_sub(u, u, qv);
        poly_free_multi(3, quo, rem, qv);
    }
    else
    {
        if (op1->deg > op2->deg)
        {
            poly *mab = poly_fast_gcd_partial_matrix(op1, op2, limit);
            poly prod1 = poly_new();
            poly prod2 = poly_new();
            poly_fast_mul(prod1, mab[2], op1);
            poly_fast_mul(prod2, mab[3], op2);
            poly_add(d, prod1, prod2);
            poly_copy(u, mab[2]);
            poly_copy(v, mab[3]);
            poly_free_multi(6, mab[0], mab[1], mab[2], mab[3], prod1, prod2);
            free(mab);
        }
        else if (op1->deg < op2->deg)
        {
            poly *mab = poly_fast_gcd_partial_matrix(op2, op1, limit);
            poly prod1 = poly_new();
            poly prod2 = poly_new();
            poly_fast_mul(prod1, mab[2], op2);
            poly_fast_mul(prod2, mab[3], op1);
            poly_add(d, prod1, prod2);
            poly_copy(u, mab[3]);
            poly_copy(v, mab[2]);
            poly_free_multi(6, mab[0], mab[1], mab[2], mab[3], prod1, prod2);
            free(mab);
        }
    }
}