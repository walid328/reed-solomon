#include <stdlib.h>
#include <stdio.h>

#include "finite_field.h"

int mod_zp(int n)
{
    int m = n % p;
    return m >= 0 ? m : m + p;
}

int add_zp(int n, int m)
{
    return mod_zp(n + m);
}

int sub_zp(int n, int m)
{
    return mod_zp(n - m);
}

int opp_zp(int n)
{
    return mod_zp(-n);
}

int mul_zp(int n, int m)
{
    return mod_zp(n * m);
}

int rand_zp()
{
    return mod_zp((rand() << 24) ^ (rand() << 16) ^ (rand() << 8) ^ rand());
}

int inv_zp(int n)
{
    if (n == 0)
    {
        fprintf(stderr, "0 is not inversible!\n");
        exit(EXIT_FAILURE);
    }
    int v0 = 0;
    int r0 = p;
    int v1 = 1;
    int r1 = n;
    while (r1 != 1)
    {
        int q = r0 / r1;
        int v2 = v0 - q * v1;
        int r2 = r0 - q * r1;
        v0 = v1;
        r0 = r1;
        v1 = v2;
        r1 = r2;
    }
    v1 %= p;
    if (v1 < 0)
        v1 += p;
    return v1;
}

int exp_zp(int base, int exp)
{
    int res = 1;
    while (exp > 0)
    {
        if (exp & 1)
            res = mul_zp(res, base);
        exp >>= 1;
        base = mul_zp(base, base);
    }
    return res;
}

int primitive_root_zp(int q, int d)
{
    int x, g;
    do
    {
        x = rand_zp();
        g = exp_zp(x, q);
    } while (x == 0 || exp_zp(g, d / 2) == 1);
    return g;
}

int min_primitive_root_zp(int q, int d)
{
	if (d == 1) return 1;
	int x = 2;
	while (exp_zp(x, d / 2) != p-1)
		x++;
    return x;
}

