#include <stdlib.h>
#include <stdio.h>

#include "finite_field.h"

int main(int argc, char **argv)
{

    int coeff_p1[] = {0, 3, 4};
    int coeff_p2[] = {0, 3, 1};

    poly *p1 = (poly *)malloc(sizeof(poly));
    p1->degree = 2;
    p1->coefficients = coeff_p1;
    poly *p2 = (poly *)malloc(sizeof(poly));
    p2->degree = 2;
    p2->coefficients = coeff_p2;
    poly *sum = add_polynomials(p1, p2);
    free(p1);
    free(p2);
    print_polynomial(sum);
    free(sum);
    return EXIT_SUCCESS;
}