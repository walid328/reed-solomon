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
#include "../src/rs_code.h"

// Tests for Reed Solomon encoding and decoding.

int p = 0;
int q = 0;
int n = 0;
zp_t omega = 0;
zp_array omegas = NULL;
zp_array inverses = NULL;

void usage(int argc, char *argv[])
{
    fprintf(stderr, "\e[1mUsage:\e[0m %s <p> <testname>\n", argv[0]);
    exit(EXIT_FAILURE);
}

void print_result(bool res, char *test_name)
{
    if (res)
        printf("Test \"%s\" finished: \e[1;32mSUCCESS\e[0m\n", test_name);
    else
        printf("Test \"%s\" finished: \e[1;31mFAILURE\e[0m\n", test_name);
}

bool test_encode(void)
{
    int block_length = 7;
    int message_length = 3;
    zp_array points = zp_array_new_set(block_length, 0, 1, 2, 3, 4, 5, 6);
    zp_array message = zp_array_new_set(message_length, 1, 2, 3);
    zp_array test_codeword = rs_encode(block_length, message_length, points, message);
    zp_array codeword = zp_array_new_set(block_length, 1, 6, 17, 34, 57, 86, 121);
    assert(zp_array_equal(test_codeword, codeword, block_length));
    zp_array_free(points);
    zp_array_free(message);
    zp_array_free(test_codeword);
    zp_array_free(codeword);
    return true;
}

bool test_fast_encode(void)
{
    int block_length = 8;
    int message_length = 3;
    zp_array message = zp_array_new_set(message_length, 1, 2, 3);
    zp_array test_codeword = rs_fast_encode(block_length, message_length, message);
    zp_array codeword = rs_encode_2(block_length, message_length, message);
    assert(zp_array_equal(test_codeword, codeword, block_length));
    zp_array_free(message);
    zp_array_free(test_codeword);
    zp_array_free(codeword);
    return true;
}

bool test_encode_decode(void)
{
    int block_length = 7;
    int message_length = 3;
    zp_array points = zp_array_new_set(block_length, 0, 1, 2, 3, 4, 5, 6);
    zp_array message = zp_array_new_set(message_length, 1, 2, 3);
    zp_array test_codeword = rs_encode(block_length, message_length, points, message);
    zp_array_add_errors(test_codeword, block_length, (block_length - message_length + 1) / 2);
    poly g_0 = poly_new();
    rs_g_0(g_0, points, block_length);
    zp_array test_message = rs_decode(g_0, block_length, message_length, points, test_codeword);
    assert(zp_array_equal(message, test_message, message_length));
    zp_array_free(points);
    zp_array_free(message);
    zp_array_free(test_codeword);
    zp_array_free(test_message);
    poly_free(g_0);
    return true;
}

bool test_encode_decode_2(void)
{
    int block_length = 8;
    int message_length = 3;
    zp_array message = zp_array_new_set(message_length, 1, 2, 3);
    zp_array test_codeword = rs_encode_2(block_length, message_length, message);
    zp_array_add_errors(test_codeword, block_length, (block_length - message_length) / 2);
    poly g_0 = poly_new();
    rs_g_0_fourier(g_0, block_length);
    zp_array test_message = rs_decode_2(g_0, block_length, message_length, test_codeword);
    assert(zp_array_equal(message, test_message, message_length));
    zp_array_free(message);
    zp_array_free(test_codeword);
    zp_array_free(test_message);
    poly_free(g_0);
    return true;
}

bool test_fast_encode_decode(void)
{
    int block_length = 8;
    int message_length = 3;
    zp_array message = zp_array_new_set(message_length, 1, 2, 3);
    zp_array test_codeword = rs_fast_encode(block_length, message_length, message);
    zp_array_add_errors(test_codeword, block_length, (block_length - message_length) / 2);
    poly g_0 = poly_new();
    rs_g_0_fourier(g_0, block_length);
    zp_array test_message = rs_fast_decode(g_0, block_length, message_length, test_codeword);
    assert(zp_array_equal(message, test_message, message_length));
    zp_array_free(message);
    zp_array_free(test_codeword);
    zp_array_free(test_message);
    poly_free(g_0);
    return true;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
        usage(argc, argv);

    srand(time(NULL));

    int p_ = atoi(argv[1]);
    field_settings_set(p_);

    bool ok;
    if (strcmp("encode", argv[2]) == 0)
        ok = test_encode();
    else if (strcmp("fast_encode", argv[2]) == 0)
        ok = test_fast_encode();
    else if (strcmp("encode_decode", argv[2]) == 0)
        ok = test_encode_decode();
    else if (strcmp("encode_decode_2", argv[2]) == 0)
        ok = test_encode_decode_2();
    else if (strcmp("fast_encode_decode", argv[2]) == 0)
        ok = test_fast_encode_decode();
    else if (strcmp("rs_code", argv[2]) == 0 || strcmp("all", argv[2]) == 0 || strcmp("*", argv[2]))
    {
        bool sub_ok = test_encode();
        print_result(sub_ok, "encode");
        ok = sub_ok;
        sub_ok = test_fast_encode();
        print_result(sub_ok, "fast_encode");
        ok &= sub_ok;
        sub_ok = test_encode_decode();
        print_result(sub_ok, "encode_decode");
        ok &= sub_ok;
        sub_ok = test_encode_decode_2();
        print_result(sub_ok, "encode_decode_2");
        ok &= sub_ok;
        sub_ok = test_fast_encode_decode();
        print_result(sub_ok, "fast_encode_decode");
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