#ifndef _XORSHIFT_H_
#define _XORSHIFT_H_

#include <stdint.h>

static uint32_t x=123456789, y=362436069, z=521288629, w=423472938;

uint32_t xorshift128(void) {
    uint32_t t = x;
    t ^= t << 11;
    t ^= t >> 8;
    x = y; y = z; z = w;
    w ^= w >> 19;
    w ^= t;
    return w;
}

#endif
