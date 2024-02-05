#include <stdlib.h>
#include <stdio.h>

#include "finite_field.h"
#include "field_settings.h"

int zp_mod(int a)
{
    int b = a % p;
    return b >= 0 ? b : b + p;
}

int zp_add(int a, int b)
{
    return zp_mod(a + b);
}

int zp_sub(int a, int b)
{
    return zp_mod(a - b);
}

int zp_opp(int a)
{
    return zp_mod(-a);
}

int zp_mul(int a, int b)
{
	long long int al = a;
	long long int bl = b;
	int ab = (al * bl) % p;
	return ab;
}

int zp_rand(void)
{
    return zp_mod((rand() << 24) ^ (rand() << 16) ^ (rand() << 8) ^ rand());
}

int zp_inv(int a)
{
    if (a == 0)
    {
        fprintf(stderr, "0 is not inversible!\n");
        exit(EXIT_FAILURE);
    }
    int v0 = 0;
    int r0 = p;
    int v1 = 1;
    int r1 = a;
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

int zp_exp(int base, int exp)
{
    int res = 1;
    while (exp > 0)
    {
        if (exp & 1)
            res = zp_mul(res, base);
        exp >>= 1;
        base = zp_mul(base, base);
    }
    return res;
}

int zp_prim_root(void)
{
    int x, g;
    do
    {
        x = zp_rand();
        g = zp_exp(x, q);
    } while (zp_exp(g, n / 2) != p - 1);
    return g;
}

int zp_prim_root_min(void)
{
    if (n == 1)
        return 1;
    int g = 2;
    while (zp_exp(g, n / 2) != p - 1)
        g++;
    return g;
}
