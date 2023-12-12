#include <stdlib.h>
#include <stdio.h>

#include "finite_field.h"

int main(int argc, char **argv)
{
    int coeff_p1[] = {5, 4, 3};
    int coeff_p2[] = {1, 3, 2};

    poly *p1 = new_poly();
    set_poly(p1, 2, coeff_p1);
    poly *p2 = new_poly();
    set_poly(p2, 2, coeff_p2);

    poly *sum = add_poly(p1, p2);
    print_poly(sum);

    poly *prod = mul_poly(p1, p2);
    print_poly(prod);

    free_poly(p1);
    free_poly(p2);
    free_poly_full(sum);
    free_poly_full(prod);
    return EXIT_SUCCESS;
}