#ifndef ASSERT_H
#define ASSERT_H

#define ASSERT(expr)                                                                        \
    do                                                                                      \
    {                                                                                       \
        if ((expr) == 0)                                                                    \
        {                                                                                   \
            fprintf(stderr, "[%s:%d] Assertion '%s' failed!\n", __FILE__, __LINE__, #expr); \
            abort();                                                                        \
        }                                                                                   \
    } while (0)

#endif