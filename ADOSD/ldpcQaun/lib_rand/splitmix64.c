/*  Written in 2015 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include "lib_rand.h"

/* This is a fixed-increment version of Java 8's SplittableRandom generator
   See http://dx.doi.org/10.1145/2714064.2660195 and
   http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

   It is a very fast generator passing BigCrush, and it can be useful if
   for some reason you absolutely want 64 bits of state; otherwise, we
   rather suggest to use a xoroshiro128+ (for moderately parallel
   computations) or xorshift1024* (for massively parallel computations)
   generator. */

static const uint64_t x0 = 0xbf3749f5b97cd3b9; /* The state can be seeded with any value. */
static const uint64_t x1 = 0x9adf8ae36c7be201;
static const uint64_t x2 = 0x8a62715ec03bedf1;
static const uint64_t x3 = 0x307e2cd5ba648f19;
static const uint64_t x4 = 0xde09275b846f1c3a;
static const uint64_t x5 = 0x0af276e98c3b514d;
static const uint64_t x6 = 0xc3a58f04d769eb21;
static const uint64_t x7 = 0xde3b278c6a5419f0;
static const uint64_t x8 = 0x9bea8164f2053c7d;
static const uint64_t x9 = 0x9b360fd2748ace51;
static const uint64_t x10 = 0x9bf382c65a7e40d1;
static const uint64_t x11 = 0x92d78f46b5ec0a13;
static const uint64_t x12 = 0x3c16fb8dea542970;
static const uint64_t x13 = 0x1a4db68ef02539c7;
static const uint64_t x14 = 0xd0964ef327c1a58b;
static const uint64_t x15 = 0x8b9f4e251c763ad0;
static const uint64_t x16 = 0x8af71e4b263c05d9;
static const uint64_t x17 = 0x0a67925eb1c84df3;
static const uint64_t x18 = 0x97d1234f8c0e6ba5;
static const uint64_t x19 = 0x1f204589ebc3d76a;
static const uint64_t x20 = 0x0ad3c86e5f9471b2;
static const uint64_t x21 = 0x5d17fc92a0683b4e;
static const uint64_t x22 = 0xb37e6a184d05c29f;
static const uint64_t x23 = 0x38247a905c1febd6;
static const uint64_t x24 = 0xfc48d925a0b1e763;
static const uint64_t x25 = 0xcadf482b36e01579;
static const uint64_t x26 = 0x2f415a376bd0ce89;
static const uint64_t x27 = 0x6458b32ca17fd90e;
static const uint64_t x28 = 0xb24a3d71cf0659e8;
static const uint64_t x29 = 0xcd30fa42876159be;
static const uint64_t x30 = 0x413fbe2c07d986a5;
static uint64_t x = x5;
uint64_t next64() {
	uint64_t z = (x += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}
