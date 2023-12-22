#include <stdlib.h>
#include <stdio.h>

#include "finite_field.h"

int add_zp(int n, int m)
{
    int sum = n + m;
    return sum < p ? sum : sum - p;
}

int sub_zp(int n, int m)
{
    int diff = n - m;
    return n >= m ? diff : diff + p;
}

int opp_zp(int n)
{
    return (p - n) % p;
}

int mul_zp(int n, int m)
{
    return (n * m) % p;
}

int rand_zp()
{
    int n = ((rand() << 24) ^ (rand() << 16) ^ (rand() << 8) ^ rand()) % p;
    if (n < 0)
        n += p;
    return n;
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