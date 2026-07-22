#pragma once
#include <algorithm>
#include <cstdint>

struct planmask {
   uint64_t w = 0 ;
   //
   // mode() = w >> 60 (high nybble, 0..15) indexes the mode tables:
   //
   //   0      : bitmask — bits 0..59 index algorithms 0..59, up to 60 activations
   //   1..5   : compact — N_SLOTS[m] sorted FIELD_W[m]-bit fields in bits 0..59,
   //            empty slot = EMPTY_V[m]; encode header = (uint64_t)m << 60
   //   6..7   : reserved (never generated)
   //   8..15  : legacy 7×9-bit — 7 sorted 9-bit fields in bits 0..62, empty=511;
   //            encode always uses m=8 (= 8ULL<<60 = 1ULL<<63 = MSB flag);
   //            resulting nybble is 8..14 or 15 depending on field-6 data
   //   ~0ULL  : none() sentinel (w==~0ULL); mode()=15 but is_none() distinguishes it
   //
   //     m  slots  bits/field  max index
   //     1    6       10        1022
   //     2    5       12        4094
   //     3    4       15       32766
   //     4    3       20     1048574
   //     5    2       30    2^30 - 2
   //   8-15    7        9         510
   //
   static constexpr int      N_SLOTS[16] = {0,6,5,4,3,2,0,0, 7,7,7,7,7,7,7,7} ;
   static constexpr int      FIELD_W[16] = {0,10,12,15,20,30,0,0, 9,9,9,9,9,9,9,9} ;
   static constexpr uint32_t EMPTY_V[16] = {0,1023,4095,32767,(1u<<20)-1,(1u<<30)-1,0,0,
                                             511,511,511,511,511,511,511,511} ;
   static constexpr uint64_t DATA_MASK   = (1ULL << 60) - 1 ;

   static planmask none()   { return { ~0ULL } ; }
   static int max_stages()  { return (1 << 30) - 1 ; }
   bool is_none()    const  { return w == ~0ULL ; }
   int  mode()       const  { return (int)(w >> 60) ; }
   bool is_bitmask() const  { return mode() == 0 ; }

   // Decode fields into out[]; return count. Works for all non-bitmask modes.
   int decode(int out[7]) const {
      int m = mode(), fw = FIELD_W[m], ns = N_SLOTS[m] ;
      uint64_t fmask = (1ULL << fw) - 1 ;
      uint32_t emp = EMPTY_V[m] ;
      int n = 0 ;
      for (int k = 0; k < ns; k++) {
         uint32_t f = (uint32_t)((w >> (k * fw)) & fmask) ;
         if (f == emp) break ;
         out[n++] = (int)f ;
      }
      return n ;
   }

   // Encode count sorted fields[] into mode m. Header = (uint64_t)m << 60;
   // for m=8 that equals 1ULL<<63, which is the legacy MSB flag.
   static planmask encode(int m, const int fields[], int count) {
      int fw = FIELD_W[m], ns = N_SLOTS[m] ;
      uint32_t emp = EMPTY_V[m] ;
      uint64_t r = (uint64_t)m << 60 ;
      for (int k = 0; k < ns; k++) {
         uint64_t v = (k < count) ? (uint32_t)fields[k] : emp ;
         r |= v << (k * fw) ;
      }
      return { r } ;
   }

   // Minimum compact mode (1-5) for index i; -1 if out of range.
   static int req_compact(int i) {
      if (i <         1023) return 1 ;
      if (i <         4095) return 2 ;
      if (i <        32767) return 3 ;
      if (i < (1<<20) - 1)  return 4 ;
      if (i < (1<<30) - 1)  return 5 ;
      return -1 ;
   }

   // Sorted insert of i into fields[0..n-1]; writes to fields[n]. Returns n+1.
   static int sorted_insert(int fields[], int n, int i) {
      int pos = n ;
      for (int k = 0; k < n; k++)
         if (i < fields[k]) { pos = k ; break ; }
      for (int k = n; k > pos; k--) fields[k] = fields[k-1] ;
      fields[pos] = i ;
      return n + 1 ;
   }

   planmask set(int i) const {
      if (is_bitmask() && i < 60) return { w | (1ULL << i) } ; // fast path (idempotent)
      if (is_none()) return none() ;
      if (test(i)) return *this ; // already set, idempotent
      int fields[61], n ; // up to 60 bitmask entries when lifting
      int cur_m = mode() ;
      if (cur_m == 0) {
         n = 0 ;
         for (int k = 0; k < 60; k++)
            if ((w >> k) & 1) fields[n++] = k ;     // lift bitmask indices
      } else {
         n = decode(fields) ;
      }
      int new_m ;
      if (i < 511) {
         // Legacy (8) preferred for small indices; stay compact if already there.
         new_m = (cur_m >= 1 && cur_m <= 5) ? cur_m : 8 ;
      } else {
         int base = (cur_m >= 1 && cur_m <= 5) ? cur_m : 0 ;
         new_m = std::max(base, req_compact(i)) ;
      }
      if (new_m < 0 || n >= N_SLOTS[new_m]) return none() ;
      return encode(new_m, fields, sorted_insert(fields, n, i)) ;
   }

   bool test(int i) const {
      if (is_bitmask()) return i < 60 && ((w >> i) & 1) ;
      if (is_none()) return false ;
      int fields[7], n = decode(fields) ;
      for (int k = 0; k < n; k++) {
         if (fields[k] == i) return true ;
         if (fields[k]  > i) return false ;
      }
      return false ;
   }

   int popcount() const {
      if (is_bitmask()) return __builtin_popcountll(w & DATA_MASK) ;
      if (is_none()) return 0 ;
      int fields[7] ;
      return decode(fields) ;
   }

   bool operator==(planmask b) const { return w == b.w ; }
   bool operator!=(planmask b) const { return w != b.w ; }
} ;
