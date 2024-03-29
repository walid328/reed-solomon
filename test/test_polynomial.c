#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

#include "../src/field_settings.h"
#include "../src/finite_field.h"
#include "../src/zp_array.h"
#include "../src/polynomial.h"

// Tests for operations on polynomials.

int p = 193;
int q = 0;
int n = 0;
zp_t omega = 0;
zp_array omegas = NULL;
zp_array inverses = NULL;

void usage(int argc, char *argv[])
{
    fprintf(stderr, "\e[1mUsage:\e[0m %s <testname>\n", argv[0]);
    exit(EXIT_FAILURE);
}

void print_result(bool res, char *test_name)
{
    if (res)
        printf("Test %s finished: \e[1;32mSUCCESS\e[0m\n", test_name);
    else
        printf("Test %s finished: \e[1;31mFAILURE\e[0m\n", test_name);
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
    poly_free_multi(2, f, g);
    return true;
}

bool test_deriv(void)
{
    poly f = poly_new_str("4*x^6 + 14*x^5 + x^3 + 13");
    poly deriv = poly_new_str("24*x^5 + 70*x^4 + 3*x^2");
    poly test_deriv = poly_new();
    poly_deriv(test_deriv, f);
    assert(poly_equal(deriv, test_deriv));
    poly_free_multi(3, f, deriv, test_deriv);
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
    poly_free_multi(4, f, g, sum, test_sum);
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
    poly_free_multi(4, f, g, prod, test_prod);
    return true;
}

bool test_mul_scalar(void)
{
    poly f = poly_new_str("3*x^5 + 3*x^4 + 7*x^3 + 4*x^2 + 4*x");
    zp_t a = 23;
    poly prod = poly_new_str("69*x^5 + 69*x^4 + 161*x^3 + 92*x^2 + 92*x");
    poly test_prod = poly_new();
    poly_mul_scalar(test_prod, a, f);
    assert(poly_equal(test_prod, prod));
    poly_free_multi(3, f, prod, test_prod);
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
    poly_free_multi(5, f, g, q, r, test_euc);
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
    poly_free_multi(8, f, g, d, u, v, uf, vg, test_xgcd);
    return true;
}

bool test_xgcd_partial(void)
{
    poly f = poly_new_rand(6);
    poly g = poly_new_rand(3);
    poly d = poly_new();
    poly u = poly_new();
    poly v = poly_new();
    int limit = 2;
    poly_xgcd_partial(d, u, v, f, g, limit);
    assert(poly_deg(d) < 2);
    poly uf = poly_new();
    poly_mul(uf, u, f);
    poly vg = poly_new();
    poly_mul(vg, v, g);
    poly test_xgcd = poly_new();
    poly_add(test_xgcd, uf, vg);
    assert(poly_equal(d, test_xgcd));
    poly_free_multi(8, f, g, d, u, v, uf, vg, test_xgcd);
    return true;
}

bool test_interpol(void)
{
    poly f = poly_new_rand(4);
    zp_array points = zp_array_new_set(5, 0, 1, 2, 3, 4);
    zp_array eval = poly_eval_zp_array(f, points, 5);
    poly test_inter = poly_new();
    interpolation(test_inter, points, eval, 5);
    assert(poly_equal(test_inter, f));
    poly_free_multi(2, f, test_inter);
    zp_array_free(points);
    zp_array_free(eval);
    return true;
}

bool test_dft(void)
{
    poly f = poly_new_rand(6);
    zp_array eval = poly_dft(f, 8);
    poly test_dft = poly_new();
    poly_inv_dft(test_dft, eval, 8);
    assert(poly_equal(test_dft, f));
    poly_free_multi(2, f, test_dft);
    zp_array_free(eval);
    return true;
}

bool test_fft(void)
{
    poly f = poly_new_rand(6);
    zp_array eval = poly_fft(f, 8);
    poly test_fft = poly_new();
    poly_inv_fft(test_fft, eval, 8);
    assert(poly_equal(test_fft, f));
    poly_free_multi(2, f, test_fft);
    zp_array_free(eval);
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
    poly_free_multi(4, f, g, prod, test_prod);
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
    poly_free_multi(5, f, g, q, r, test_euc);
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
    poly_free_multi(8, f, g, d, u, v, uf, vg, test_xgcd);
    return true;
}

bool test_fast_xgcd_partial(void)
{
    poly f = poly_new_rand(6);
    poly g = poly_new_rand(3);
    poly d = poly_new();
    poly u = poly_new();
    poly v = poly_new();
    int limit = 2;
    poly_fast_xgcd_partial(d, u, v, f, g, limit);
    assert(poly_deg(d) < 2);
    poly uf = poly_new();
    poly_mul(uf, u, f);
    poly vg = poly_new();
    poly_mul(vg, v, g);
    poly test_xgcd = poly_new();
    poly_add(test_xgcd, uf, vg);
    assert(poly_equal(d, test_xgcd));
    poly_free_multi(8, f, g, d, u, v, uf, vg, test_xgcd);
    return true;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
        usage(argc, argv);

    srand(time(NULL));

    int p_ = atoi(argv[1]);
    field_settings_set(p_);

    // start test
    bool ok;
    if (strcmp("new_free", argv[2]) == 0)
        ok = test_new_free();
    else if (strcmp("add", argv[2]) == 0)
        ok = test_add();
    else if (strcmp("mul", argv[2]) == 0)
        ok = test_mul();
    else if (strcmp("copy", argv[2]) == 0)
        ok = test_copy();
    else if (strcmp("deriv", argv[2]) == 0)
        ok = test_deriv();
    else if (strcmp("mul_scalar", argv[2]) == 0)
        ok = test_mul_scalar();
    else if (strcmp("euc_div", argv[2]) == 0)
        ok = test_euc_div();
    else if (strcmp("xgcd", argv[2]) == 0)
        ok = test_xgcd();
    else if (strcmp("xgcd_partial", argv[2]) == 0)
        ok = test_xgcd_partial();
    else if (strcmp("interpol", argv[2]) == 0)
        ok = test_interpol();
    else if (strcmp("dft", argv[2]) == 0)
        ok = test_dft();
    else if (strcmp("fft", argv[2]) == 0)
        ok = test_fft();
    else if (strcmp("fast_mul", argv[2]) == 0)
        ok = test_fast_mul();
    else if (strcmp("fast_euc_div", argv[2]) == 0)
        ok = test_fast_euc_div();
    else if (strcmp("fast_xgcd", argv[2]) == 0)
        ok = test_fast_xgcd();
    else if (strcmp("fast_xgcd_partial", argv[2]) == 0)
        ok = test_fast_xgcd_partial();
    else if (strcmp("polynomial", argv[2]) == 0 || strcmp("all", argv[2]) == 0 || strcmp("*", argv[2]) == 0)
    {
        bool sub_ok = test_new_free();
        print_result(sub_ok, "new_free");
        ok = sub_ok;
        sub_ok = test_add();
        print_result(sub_ok, "add");
        ok &= sub_ok;
        sub_ok = test_mul();
        print_result(sub_ok, "mul");
        ok &= sub_ok;
        sub_ok = test_copy();
        print_result(sub_ok, "copy");
        ok &= sub_ok;
        sub_ok = test_deriv();
        print_result(sub_ok, "deriv");
        ok &= sub_ok;
        sub_ok = test_mul_scalar();
        print_result(sub_ok, "mul_scalar");
        ok &= sub_ok;
        sub_ok = test_euc_div();
        print_result(sub_ok, "euc_div");
        ok &= sub_ok;
        sub_ok = test_xgcd();
        print_result(sub_ok, "xgcd");
        ok &= sub_ok;
        sub_ok = test_xgcd_partial();
        print_result(sub_ok, "xgcd_partial");
        ok &= sub_ok;
        sub_ok = test_interpol();
        print_result(sub_ok, "interpol");
        ok &= sub_ok;
        sub_ok = test_dft();
        print_result(sub_ok, "dft");
        ok &= sub_ok;
        sub_ok = test_fft();
        print_result(sub_ok, "fft");
        ok &= sub_ok;
        sub_ok = test_fast_mul();
        print_result(sub_ok, "fast_mul");
        ok &= sub_ok;
        sub_ok = test_fast_euc_div();
        print_result(sub_ok, "fast_euc_div");
        ok &= sub_ok;
        sub_ok = test_fast_xgcd();
        print_result(sub_ok, "fast_xgcd");
        ok &= sub_ok;
        sub_ok = test_fast_xgcd_partial();
        print_result(sub_ok, "fast_xgcd_partial");
        ok &= sub_ok;
    }

    else
    {
        fprintf(stderr, "Error: test \"%s\" not found!\n", argv[2]);
        exit(EXIT_FAILURE);
    }

    printf("==> ");
    print_result(ok, argv[2]);

    field_settings_free();
    return EXIT_SUCCESS;
}