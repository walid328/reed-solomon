#include <stdlib.h>
#include <stdio.h>

#include "finite_field.h"
#include "field_settings.h"

zp_t zp_add(zp_t a, zp_t b)
{
    zp_t sum = a + b;
    if (sum >= p)
        return sum - p;
    return sum;
}

zp_t zp_sub(zp_t a, zp_t b)
{
    if (b > a)
        return (a + p) - b;
    return a - b;
}

zp_t zp_opp(zp_t a)
{
    if (a == 0)
        return 0;
    return p - a;
}

zp_t zp_mul(zp_t a, zp_t b)
{
    unsigned long long int al = a;
    unsigned long long int bl = b;
    zp_t ab = (al * bl) % p;
    return ab;
}

zp_t zp_rand(void)
{
    zp_t r = (rand() << 24) ^ (rand() << 16) ^ (rand() << 8) ^ rand();
    return r % p;
}

zp_t zp_inv_xgcd(zp_t a)
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

zp_t zp_inv(zp_t a)
{
    return inverses[a];
}

zp_t zp_exp(zp_t base, zp_t exp)
{
    zp_t res = 1;
    while (exp > 0)
    {
        if (exp & 1)
            res = zp_mul(res, base);
        exp >>= 1;
        base = zp_mul(base, base);
    }
    return res;
}

zp_t zp_prim_root(void)
{
    zp_t x, g;
    do
    {
        x = zp_rand();
        g = zp_exp(x, q);
    } while (zp_exp(g, n / 2) != p - 1);
    return g;
}

zp_t zp_prim_root_min(void)
{
    if (n == 1)
        return 1;
    zp_t g = 2;
    while (zp_exp(g, n / 2) != p - 1)
        g++;
    return g;
}
