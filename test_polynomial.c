#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "field_settings.h"
#include "finite_field.h"
#include "array.h"
#include "polynomial.h"

int p = 193;
int q = 0;
int n = 0;
int omega = 0;
array omegas = NULL;

void usage(int argc, char *argv[])
{
    fprintf(stderr, "Usage: %s <testname> [<...>]\n", argv[0]);
    exit(EXIT_FAILURE);
}

bool test_new_free(void)
{
    poly f = poly_new();
    assert(f);
    poly_free(f);
    return true;
}

bool test_copy(void)
{
    poly f = poly_new_rand(10);
    poly g = poly_new();
    poly_copy(g, f);
    assert(poly_equal(f, g));
    poly_free(f);
    poly_free(g);
    return true;
}

bool test_rev(void)
{
    poly f = poly_new_str("4*x^4 + 2*x^2");
    poly rev = poly_new_str("2*x^2 + 4");
    poly test_rev = poly_new();
    poly_rev(test_rev, f);
    assert(poly_equal(rev, test_rev));
    poly_free(f);
    poly_free(rev);
    poly_free(test_rev);
    return true;
}

bool test_deriv(void)
{
    poly f = poly_new_str("4*x^6 + 14*x^5 + x^3 + 13");
    poly deriv = poly_new_str("24*x^5 + 70*x^4 + 3*x^2");
    poly test_deriv = poly_new();
    poly_deriv(test_deriv, f);
    assert(poly_equal(deriv, test_deriv));
    poly_free(f);
    poly_free(deriv);
    poly_free(test_deriv);
    return true;
}

bool test_add(void)
{
    poly f = poly_new_str("x^2 + x + 1");
    poly g = poly_new_str("3*x^3 + 4*x");
    poly sum = poly_new_str("3*x^3 + x^2 + 5*x + 1");
    poly test_sum = poly_new();
    poly_add(test_sum, f, g);
    assert(poly_equal(test_sum, sum));
    poly_free(f);
    poly_free(g);
    poly_free(sum);
    poly_free(test_sum);
    return true;
}

bool test_mul(void)
{
    poly f = poly_new_str("x^2 + x + 1");
    poly g = poly_new_str("3*x^3 + 4*x");
    poly prod = poly_new_str("3*x^5 + 3*x^4 + 7*x^3 + 4*x^2 + 4*x");
    poly test_prod = poly_new();
    poly_mul(test_prod, f, g);
    assert(poly_equal(test_prod, prod));
    poly_free(f);
    poly_free(g);
    poly_free(prod);
    poly_free(test_prod);
    return true;
}

bool test_mul_scalar(void)
{
    poly f = poly_new_str("3*x^5 + 3*x^4 + 7*x^3 + 4*x^2 + 4*x");
    int a = 23;
    poly prod = poly_new_str("69*x^5 + 69*x^4 + 161*x^3 + 92*x^2 + 92*x");
    poly test_prod = poly_new();
    poly_mul_scalar(test_prod, a, f);
    assert(poly_equal(test_prod, prod));
    poly_free(f);
    poly_free(prod);
    poly_free(test_prod);
    return true;
}

bool test_euc_div(void)
{
    poly f = poly_new_rand(6);
    poly g = poly_new_rand(3);
    poly q = poly_new();
    poly r = poly_new();
    poly_euc_div(q, r, f, g);
    assert(poly_deg(r) < poly_deg(g));
    poly test_euc = poly_new();
    poly_mul(test_euc, g, q);
    poly_add(test_euc, test_euc, r);
    assert(poly_equal(test_euc, f));
    poly_free(f);
    poly_free(g);
    poly_free(q);
    poly_free(r);
    poly_free(test_euc);
    return true;
}

bool test_xgcd(void)
{
    poly f = poly_new_rand(6);
    poly g = poly_new_rand(3);
    poly d = poly_new();
    poly u = poly_new();
    poly v = poly_new();
    poly_xgcd(d, u, v, f, g);
    poly uf = poly_new();
    poly_mul(uf, u, f);
    poly vg = poly_new();
    poly_mul(vg, v, g);
    poly test_xgcd = poly_new();
    poly_add(test_xgcd, uf, vg);
    assert(poly_equal(d, test_xgcd));
    poly_free(f);
    poly_free(g);
    poly_free(d);
    poly_free(u);
    poly_free(v);
    poly_free(uf);
    poly_free(vg);
    poly_free(test_xgcd);
    return true;
}

bool test_interpol(void)
{
    poly f = poly_new_rand(4);
    array points = array_new_set(5, 0, 1, 2, 3, 4);
    array eval = poly_eval_array(f, points, 5);
    poly test_inter = poly_new();
    interpolation(test_inter, points, eval, 5);
    assert(poly_equal(test_inter, f));
    poly_free(f);
    poly_free(test_inter);
    array_free(points);
    array_free(eval);
    return true;
}

bool test_dft(void)
{
    poly f = poly_new_rand(6);
    array eval = poly_dft(f);
    poly test_dft = poly_new();
    poly_inv_dft(test_dft, eval);
    assert(poly_equal(test_dft, f));
    poly_free(f);
    poly_free(test_dft);
    array_free(eval);
    return true;
}

bool test_fft(void)
{
    poly f = poly_new_rand(6);
    array eval = poly_fft(f);
    poly test_fft = poly_new();
    poly_inv_fft(test_fft, eval);
    assert(poly_equal(test_fft, f));
    poly_free(f);
    poly_free(test_fft);
    array_free(eval);
    return true;
}

bool test_fast_mul(void)
{
    poly f = poly_new_str("x^2 + x + 1");
    poly g = poly_new_str("3*x^3 + 4*x");
    poly prod = poly_new_str("3*x^5 + 3*x^4 + 7*x^3 + 4*x^2 + 4*x");
    poly test_prod = poly_new();
    poly_fast_mul(test_prod, f, g);
    assert(poly_equal(test_prod, prod));
    poly_free(f);
    poly_free(g);
    poly_free(prod);
    poly_free(test_prod);
    return true;
}

bool test_fast_euc_div(void)
{
    poly f = poly_new_rand(6);
    poly g = poly_new_rand(3);
    poly q = poly_new();
    poly r = poly_new();
    poly_fast_euc_div(q, r, f, g);
    assert(poly_deg(r) < poly_deg(g));
    poly test_euc = poly_new();
    poly_mul(test_euc, g, q);
    poly_add(test_euc, test_euc, r);
    assert(poly_equal(test_euc, f));
    poly_free(f);
    poly_free(g);
    poly_free(q);
    poly_free(r);
    poly_free(test_euc);
    return true;
}

bool test_fast_xgcd(void)
{
    poly f = poly_new_rand(6);
    poly g = poly_new_rand(3);
    poly d = poly_new();
    poly u = poly_new();
    poly v = poly_new();
    poly_fast_xgcd(d, u, v, f, g);
    poly uf = poly_new();
    poly_mul(uf, u, f);
    poly vg = poly_new();
    poly_mul(vg, v, g);
    poly test_xgcd = poly_new();
    poly_add(test_xgcd, uf, vg);
    assert(poly_equal(d, test_xgcd));
    poly_free(f);
    poly_free(g);
    poly_free(d);
    poly_free(u);
    poly_free(v);
    poly_free(uf);
    poly_free(vg);
    poly_free(test_xgcd);
    return true;
}

int main(int argc, char *argv[])
{
    field_settings_update();

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
    else if (strcmp("copy", argv[1]) == 0)
        ok = test_copy();
    else if (strcmp("rev", argv[1]) == 0)
        ok = test_rev();
    else if (strcmp("deriv", argv[1]) == 0)
        ok = test_deriv();
    else if (strcmp("mul_scalar", argv[1]) == 0)
        ok = test_mul_scalar();
    else if (strcmp("euc_div", argv[1]) == 0)
        ok = test_euc_div();
    else if (strcmp("xgcd", argv[1]) == 0)
        ok = test_xgcd();
    else if (strcmp("interpol", argv[1]) == 0)
        ok = test_interpol();
    else if (strcmp("dft", argv[1]) == 0)
        ok = test_dft();
    else if (strcmp("fft", argv[1]) == 0)
        ok = test_fft();
    else if (strcmp("fast_mul", argv[1]) == 0)
        ok = test_fast_mul();
    else if (strcmp("fast_euc_div", argv[1]) == 0)
        ok = test_fast_euc_div();
    else if (strcmp("fast_xgcd", argv[1]) == 0)
        ok = test_fast_xgcd();

    else
    {
        fprintf(stderr, "Error: test \"%s\" not found!\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    // print test result
    if (ok)
    {
        fprintf(stderr, "Test \"%s\" finished: \033[0;32mSUCCESS\033[0m\n", argv[1]);
        field_settings_free();
        return EXIT_SUCCESS;
    }
    else
    {
        fprintf(stderr, "Test \"%s\" finished: \033[0;31mFAILURE\033[0m\n", argv[1]);
        field_settings_free();
        return EXIT_FAILURE;
    }
}