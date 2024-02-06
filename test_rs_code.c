#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "field_settings.h"
#include "finite_field.h"
#include "array.h"
#include "polynomial.h"
#include "rs_code.h"

int p = 0;
int q = 0;
int n = 0;
int omega = 0;
array omegas = NULL;
array inverses = NULL;

void usage(int argc, char *argv[])
{
    fprintf(stderr, "Usage: %s <testname> [<...>]\n", argv[0]);
    exit(EXIT_FAILURE);
}

bool test_encode(void)
{
    int block_length = 7;
    int message_length = 3;
    array points = array_new_set(block_length, 0, 1, 2, 3, 4, 5, 6);
    array message = array_new_set(message_length, 1, 2, 3);
    array test_codeword = rs_encode(block_length, message_length, points, message);
    array codeword = array_new_set(block_length, 1, 6, 17, 34, 57, 86, 121);
    assert(array_equal(test_codeword, codeword, block_length));
    array_free(points);
    array_free(message);
    array_free(test_codeword);
    array_free(codeword);
    return true;
}

bool test_fast_encode(void)
{
    int block_length = 8;
    int message_length = 3;
    array message = array_new_set(message_length, 1, 2, 3);
    array test_codeword = rs_fast_encode(block_length, message_length, message);
    array codeword = rs_encode_2(block_length, message_length, message);
    assert(array_equal(test_codeword, codeword, block_length));
    array_free(message);
    array_free(test_codeword);
    array_free(codeword);
    return true;
}

bool test_decode(void)
{
    int block_length = 7;
    int message_length = 3;
    array points = array_new_set(block_length, 0, 1, 2, 3, 4, 5, 6);
    array received = array_new_set(block_length, 1, 6, 123, 456, 57, 86, 121);
    array test_message = rs_decode(block_length, message_length, points, received);
    array message = array_new_set(message_length, 1, 2, 3);
    assert(array_equal(test_message, message, message_length));
    array_free(points);
    array_free(received);
    array_free(test_message);
    array_free(message);
    return true;
}

bool test_encode_decode(void)
{
    int block_length = 7;
    int message_length = 3;
    array points = array_new_set(n, 0, 1, 2, 3, 4, 5, 6);
    array message = array_new_set(message_length, 1, 2, 3);
    array test_codeword = rs_encode(block_length, message_length, points, message);
    array_add_errors(test_codeword, block_length, (block_length - message_length + 1) / 2);
	array test_message = rs_decode(block_length, message_length, points, test_codeword);
    assert(array_equal(message, test_message, message_length));
    array_free(points);
    array_free(message);
    array_free(test_codeword);
    array_free(test_message);
    return true;
}

bool test_encode_decode_2(void)
{
    int block_length = 8;
    int message_length = 3;
    array message = array_new_set(message_length, 1, 2, 3);
    array test_codeword = rs_encode_2(block_length, message_length, message);
    array_add_errors(test_codeword, block_length, (block_length - message_length + 1) / 2);
    array test_message = rs_decode_2(block_length, message_length, test_codeword);
    assert(array_equal(message, test_message, message_length));
    array_free(message);
    array_free(test_codeword);
    array_free(test_message);
    return true;
}

bool test_fast_encode_decode(void)
{
    int block_length = 8;
    int message_length = 3;
    array message = array_new_set(message_length, 1, 2, 3);
    array test_codeword = rs_fast_encode(block_length, message_length, message);
    array_add_errors(test_codeword, block_length, (block_length - message_length + 1) / 2);
    array test_message = rs_fast_decode(block_length, message_length, test_codeword);
    assert(array_equal(message, test_message, message_length));
    array_free(message);
    array_free(test_codeword);
    array_free(test_message);
    return true;
}

int main(int argc, char *argv[])
{
    field_settings_set(193);

    if (argc == 1)
        usage(argc, argv);

    // start test
    fprintf(stderr, "=> Start test \"%s\"\n", argv[1]);
    bool ok = false;
    if (strcmp("encode", argv[1]) == 0)
        ok = test_encode();
    else if (strcmp("fast_encode", argv[1]) == 0)
        ok = test_fast_encode();
    else if (strcmp("decode", argv[1]) == 0)
        ok = test_decode();
    else if (strcmp("encode_decode", argv[1]) == 0)
        ok = test_encode_decode();
    else if (strcmp("encode_decode_2", argv[1]) == 0)
        ok = test_encode_decode_2();
    else if (strcmp("fast_encode_decode", argv[1]) == 0)
        ok = test_fast_encode_decode();

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