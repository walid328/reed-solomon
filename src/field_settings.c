#include <stdlib.h>
#include <assert.h>

#include "field_settings.h"

void field_settings_set(int p_)
{
    field_settings_free();
    p = p_;
    q = p - 1;
    n = 1;
    while ((q & 1) == 0)
    {
        q >>= 1;
        n <<= 1;
    }
    omega = zp_prim_root_min();
    omegas = zp_array_new(n);
    omegas[0] = 1;
    for (int i = 1; i < n; i++)
    {
        omegas[i] = zp_mul(omegas[i - 1], omega);
    }
    inverses = zp_array_new(p);
    inverses[0] = 0;
    for (int i = 1; i < p; i++)
        inverses[i] = zp_inv_xgcd(i);
}

void field_settings_free(void)
{
    zp_array_free(omegas);
    zp_array_free(inverses);
}