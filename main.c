#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "src/field_settings.h"
#include "src/finite_field.h"
#include "src/zp_array.h"
#include "src/polynomial.h"
#include "src/rs_code.h"

// Encode and decode using Reed Solomon codes from the command line.

int p = 0;
int q = 0;
int n = 0;
zp_t omega = 0;
zp_array omegas = NULL;
zp_array inverses = NULL;

void usage(int argc, char *argv[])
{
    fprintf(stderr, "\e[1;31mFailure: \e[0mNot enough parameters!\n");
    fprintf(stderr, "\e[1mUsage:\e[0m %s -option [args]\n", argv[0]);
    fprintf(stderr, "Use -h for more information.\n");
    exit(EXIT_FAILURE);
}

// Print in sdtout the codeword associate to a given message,
// using RS(n,k) with fast operations.
void encode(int argc, char *argv[])
{
    int block_length = atoi(argv[3]);
    int message_length = atoi(argv[4]);
    zp_array message = zp_array_new(message_length);
    for (int i = 0; i < message_length; i++)
        message[i] = atoi(argv[5 + i]);
    zp_array codeword = rs_fast_encode(block_length, message_length, message);
    printf("The codeword is:\n");
    zp_array_print(codeword, block_length);
    zp_array_free(message);
    zp_array_free(codeword);
}

// Print in sdtout the message associate to a given received,
// using RS(n,k) with fast operations.
void decode(int argc, char *argv[])
{
    int block_length = atoi(argv[3]);
    int message_length = atoi(argv[4]);
    zp_array received = zp_array_new(block_length);
    for (int i = 0; i < block_length; i++)
        received[i] = atoi(argv[5 + i]);
    poly g_0 = poly_new();
    rs_g_0_fourier(g_0, block_length);
    zp_array message = rs_fast_decode(g_0, block_length, message_length, received);
    printf("The message is:\n");
    zp_array_print(message, message_length);
    zp_array_free(received);
    zp_array_free(message);
    poly_free(g_0);
}

int main(int argc, char *argv[])
{
    if (argc == 1)
        usage(argc, argv);
    else if (strcmp("-h", argv[1]) == 0)
    {
        printf("options\n");
        printf("-h \t\t\t\t show this message\n");
        printf("-e <p> <n> <k> <message>  \t encoding messsage using RS(n,k) over GF(p)\n");
        printf("-d <p> <n> <k> <received> \t decoding received using RS(n,k) over GF(p)\n");
        printf("==> message and received should be represented as a list with elements separated by a space.\n");
        printf("==> user should be careful to enter the right amount of elements.\n");
        printf("==> message should have k elements, received should have n elements.\n");
    }
    else
    {
        if (argc < 5)
            usage(argc, argv);
        int p_ = atoi(argv[2]);
        field_settings_set(p_);
        if (strcmp("-e", argv[1]) == 0)
            encode(argc, argv);
        else if (strcmp("-d", argv[1]) == 0)
            decode(argc, argv);
        field_settings_free();
    }
    return EXIT_SUCCESS;
}