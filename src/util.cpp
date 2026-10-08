#include "util.h"
#include <iostream>
#include <cstdlib>

std::mt19937_64 rng(std::random_device{}()) ;

static uint64_t splitmix64(uint64_t x) {
   x += 0x9e3779b97f4a7c15ULL ;
   x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL ;
   x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL ;
   return x ^ (x >> 31) ;
}
uint64_t mixseed(uint64_t a, uint64_t b, uint64_t c) {
   return splitmix64(splitmix64(splitmix64(a) ^ b) ^ c) ;
}

void error(const char *s) {
   std::cerr << s << std::endl ;
   std::exit(10) ;
}
