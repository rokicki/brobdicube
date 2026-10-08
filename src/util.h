#pragma once
#include <random>

extern std::mt19937_64 rng ;
// splitmix64-style mix of three values into one well-distributed seed.
uint64_t mixseed(uint64_t a, uint64_t b, uint64_t c) ;
void error(const char *s) ;
