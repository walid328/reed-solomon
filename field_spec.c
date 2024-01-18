#include <stdlib.h>

#include "field_spec.h"
#include "finite_field.h"

void update_field_spec(void)
{
    q = p - 1;
    n = 1;
    while ((q & 1) == 0)
    {
        q >>= 1;
        n <<= 1;
    }
    omega = zp_prim_root_min();
    omegas = (int *)malloc((n + 1) * sizeof(int));
    omegas[0] = 1;
    for (int i = 1; i <= n; i++)
    {
        omegas[i] = zp_mul(omegas[i - 1], omega);
    }
}