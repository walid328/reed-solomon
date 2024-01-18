#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "field_spec.h"
#include "polynomial.h"
#include "finite_field.h"

int p = 193;
int q = 0;
int n = 0;
int omega = 0;
int *omegas = NULL;

#define ASSERT(expr)                                                                        \
    do                                                                                      \
    {                                                                                       \
        if ((expr) == 0)                                                                    \
        {                                                                                   \
            fprintf(stderr, "[%s:%d] Assertion '%s' failed!\n", __FILE__, __LINE__, #expr); \
            abort();                                                                        \
        }                                                                                   \
    } while (0)

void usage(int argc, char *argv[])
{
    fprintf(stderr, "Usage: %s <testname> [<...>]\n", argv[0]);
    exit(EXIT_FAILURE);
}

bool poly_cmp(poly *f, poly *g)
{
    if (poly_degree(f) != poly_degree(g))
        return false;
    for (int i = 0; i < poly_degree(f); i++)
        if (poly_coeffs(f)[i] != poly_coeffs(g)[i])
            return false;
    return true;
}

bool test_new_free(void)
{
    poly *f = poly_new();
    ASSERT(f);
    poly_free(f);
    return true;
}

bool test_add(void)
{
    poly *f = poly_new_from_str("x^2 + x + 1");
    ASSERT(f);
    poly *g = poly_new_from_str("3*x^3 + 4*x");
    ASSERT(g);
    poly *test_sum = poly_new();
    ASSERT(test_sum);
    poly_add(test_sum, f, g);
    poly *sum = poly_new_from_str("3*x^3 + x^2 + 5*x + 1");
    ASSERT(sum);
    ASSERT(poly_cmp(test_sum, sum));
    poly_free_full(f);
    poly_free_full(g);
    poly_free_full(sum);
    poly_free_full(test_sum);
    return true;
}

bool test_mul(void)
{
    poly *f = poly_new_from_str("x^2 + x + 1");
    ASSERT(f);
    poly *g = poly_new_from_str("3*x^3 + 4*x");
    ASSERT(g);
    poly *test_prod = poly_new();
    ASSERT(test_prod);
    poly_mul(test_prod, f, g);
    poly *prod = poly_new_from_str("3*x^5 + 3*x^4 + 7*x^3 + 4*x^2 + 4*x");
    ASSERT(prod);
    ASSERT(poly_cmp(test_prod, prod));
    poly_free_full(f);
    poly_free_full(g);
    poly_free_full(prod);
    poly_free_full(test_prod);
    return true;
}

int main(int argc, char *argv[])
{
    update_field_spec();

    if (argc == 1)
        usage(argc, argv);

    // start test
    fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
    bool ok = false;
    if (strcmp("new_free", argv[1]) == 0)
        ok = test_new_free();
    else if (strcmp("add", argv[1]) == 0)
        ok = test_add();
    else if (strcmp("mul", argv[1]) == 0)
        ok = test_mul();

    else
    {
        fprintf(stderr, "Error: test \"%s\" not found!\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    // print test result
    if (ok)
    {
        fprintf(stderr, "Test \"%s\" finished: SUCCESS\n", argv[1]);
        free(omegas);
        return EXIT_SUCCESS;
    }
    else
    {
        fprintf(stderr, "Test \"%s\" finished: FAILURE\n", argv[1]);
        free(omegas);
        return EXIT_FAILURE;
    }
}