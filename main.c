
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "field_settings.h"
#include "finite_field.h"
#include "array.h"
#include "polynomial.h"
#include "rs_code.h"

// Encode and decode using Reed Solomon codes from the command line.

int p = 0;
int q = 0;
int n = 0;
int omega = 0;
array omegas = NULL;

void usage(int argc, char *argv[])
{
    fprintf(stderr, "\e[1;31mFailure: \e[0mNot enough parameters!\nUsage: %s -option [args]\nUse -h for more information.\n", argv[0]);
    exit(EXIT_FAILURE);
}

// Print in sdtout the codeword associate to a given message,
// using RS(n,k) with fast operations.
void encode(int argc, char *argv[])
{
    int block_length = atoi(argv[3]);
    int message_length = atoi(argv[4]);
    array message = array_new(message_length);
    for (int i = 0; i < message_length; i++)
        message[i] = atoi(argv[5 + i]);
    array codeword = rs_fast_encode(block_length, message_length, message);
    printf("The codeword is:\n");
    array_print(codeword, block_length);
    array_free(message);
    array_free(codeword);
}

// Print in sdtout the message associate to a given received,
// using RS(n,k) with fast operations.
void decode(int argc, char *argv[])
{
    int block_length = atoi(argv[3]);
    int message_length = atoi(argv[4]);
    array received = array_new(block_length);
    for (int i = 0; i < block_length; i++)
        received[i] = atoi(argv[5 + i]);
    array message = rs_fast_decode(block_length, message_length, received);
    printf("The message is:\n");
    array_print(message, message_length);
    array_free(received);
    array_free(message);
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