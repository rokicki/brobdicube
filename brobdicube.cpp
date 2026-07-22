#ifndef WRITEFRAMES
#define WRITEFRAMES (0)
#endif
#include <algorithm>
#include <array>
#include <assert.h>
#include <strings.h>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <random>
#include <map>
#include <queue>
#include <vector>
#include <dirent.h>
#include <sys/time.h>
#include <thread>
#include <mutex>
#include <atomic>
#include <memory>
#if WRITEFRAMES
#include "fpng.h"
#endif
using namespace std ;
using ll = long long ;
using lll = __int128_t ;
// In the following like, you can change ll to lll so you can run
// estimates up to about 10^34 or so (until something overflows).
// But this may hurt performance, so by default we leave it at ll.
using cubesizetype = ll ;
static cubesizetype atocsz(const char *s) {
   cubesizetype r = 0 ;
   int neg = (*s == '-') ; if (neg) s++ ;
   while (*s >= '0' && *s <= '9') r = r * 10 + (*s++ - '0') ;
   return neg ? -r : r ;
}
#include "planmask.h"
#include "xorbit.h"
#include "mv.h"
int showdetails = 0 ;
int alwaysinmem = 0 ;
ll writemoveinc ;
ll writemovetarget = 0 ;
ll runningmovecount = 0 ;
int fast = 1 ;
int readmoves, writemoves ;
ll randmoves = 0 ;
ll cancellations = 0 ;
ll blocks = 0 ;
ll wingmoves = 0 ; 
int extrapolate = 0 ;
int stats = 0 ;
double extrapolatemul = 1 ;
ll solvecount = 1 ;
cubesizetype N ;
cubesizetype phasewidth = -1 ;
int planonly = 0 ;
int useavg = 1 ;
int greedyonly = 0 ;
cubesizetype plan_range_start = 0, plan_range_end = 0 ;
string algsdir = "./algs3" ;
int wing_cycle_index[24][24][24] ;   // -1 if the triple is not a valid 3-cycle, inverse of wing_cycle_triple
int wing_cycle_triple[4048][3] ;     // inverse of wing_cycle_index 
const int oppface[] = { 5, 3, 4, 1, 2, 0 } ;
const char *faceorder = "ULFRBD" ;
#include "util.h"
static double start ;
static mutex io_mutex ;
static int nprob_threads = (int)thread::hardware_concurrency() ;
static bool describealgorithms = false ;
double walltime() {
   struct timeval tv ;
   gettimeofday(&tv, 0) ;
   return tv.tv_sec + 0.000001 * tv.tv_usec ;
}
double duration() {
   double now = walltime() ;
   double r = now - start ;
   start = now ;
   return r ;
}
ostream& operator<<(ostream& os, __int128_t n) {
    if (n == 0) return os << "0";
    if (n < 0) {
        os << "-";
        n = -n;
    }
    string s;
    while (n > 0) {
        s += (char)('0' + (n % 10));
        n /= 10;
    }
    reverse(s.begin(), s.end());
    return os << s;
}
/*
 *   These are the 24 rotations, given as face mappings.  We only give the
 *   first three and build the rest.
 */
int rots[24][6] = {
   {0, 1, 2, 3, 4, 5}, {0, 2, 3, 4, 1, 5}, {2, 1, 5, 3, 0, 4},
} ;
void coordtospeffz(int Mx, int My, int x, int y, int &ox, int &oy, int &oz) {
   assert(x >= 1) ;
   assert(x <= Mx + My) ;
   assert(y >= 1) ;
   assert(y <= Mx + My) ;
   if (y <= My && x <= Mx) {
      ox = x ;
      oy = y ;
      oz = 0 ;
   } else if (y <= Mx && x > Mx) {
      oy = Mx + My + 1 - x ;
      ox = y ;
      oz = 1 ;
   } else if (y > Mx && x > My) {
      ox = Mx + My + 1 - x ;
      oy = Mx + My + 1 - y ;
      oz = 2 ;
   } else if (y > My && x <= My) {
      ox = Mx + My + 1 - y ;
      oy = x ;
      oz = 3 ;
   } else {
      ox = 1 ;
      oy = 1 ;
      oz = 4 ; // special indicator for central faces
   }
   assert(ox <= Mx) ;
   assert(oy <= My) ;
   assert(ox >= 1) ;
   assert(oy >= 1) ;
}
vector<mv> moveset ;
struct xcube ;
void domove(mv mv, xcube &cube, bool suppress_cancellations=false) ;
#include "threecube.h"
template<int> int goodmove(int, int) { return 1 ; }
template<> int goodmove<1>(int f, int tw) {
   return f==0 || f==5 || tw == 2 ;
}
template<int dep, int size, int maxdepth, int limitmoves=0> struct prunetable {
   prunetable(const char *name_) : solved('C', 0), name(name_) {
      memset(prune, dep, size) ;
   }
   void filltable(threecube &cc, int togo, int prevf, int d) {
      if (togo == 0) {
         int h = cc.hash() % size ;
         if (prune[h] == dep) {
            prune[h] = d ;
         }
      } else {
         for (int f=0; f<6; f++)
            if (f != prevf) {
               for (int tw=1; tw<4; tw++) {
                  cc.domove({1, (short)f, 1}) ;
                  if (goodmove<limitmoves>(f, tw) &&
                      (d != togo || cc != solved))
                     filltable(cc, togo-1, f, d) ;
               }
               cc.domove({1, (short)f, 1}) ;
            }
      }
   }
   void fillit(threecube &cc) {
      duration() ;
      solved = cc ;
      for (int d=0; d<dep; d++)
         filltable(cc, d, -1, d) ;
      cout << "Created pruning table " << name << " in " << duration() << endl ;
   }
   int lookup(const threecube &cc) {
      return prune[cc.hash() % size] ;
   }
   int solve(threecube &cc, int togo, int prevf) {
      if (togo == 0)
         return (cc == solved) ;
      if (lookup(cc) > togo)
         return 0 ;
      for (int f=0; f<6; f++)
         if (f != prevf) {
            for (int tw=1; tw<4; tw++) {
               mvs.push_back({1, (short)f, (short)tw}) ;
               cc.domove({1, (short)f, 1}) ;
               if (goodmove<limitmoves>(f, tw) && solve(cc, togo-1, f))
                  return 1 ;
               mvs.pop_back() ;
            }
            cc.domove({1, (short)f, 1}) ;
         }
      return 0 ;
   }
   int solve(threecube &cc) {
      duration() ;
      mvs.clear() ;
      for (int curdep=lookup(cc); curdep<=maxdepth; curdep++) {
         if (solve(cc, curdep, -1)) {
            cout << "Solved " << name << " in " << duration() << endl ;
            return curdep ;
         }
      }
      return -1 ;
   }
   unsigned char prune[size] ;
   threecube solved ;
   vector<mv> mvs ;
   const char *name ;
} ;
struct xcube {
   xcube(cubesizetype Mx_, cubesizetype My_, int src_) :
  Mx(Mx_), My(My_), rows(My), tc(src_, (Mx_ + My_) & 1), src(src_) {
      for (int i=0; i<6; i++)
         facerots[i] = 0 ;
      // central faces, only for odd cubes
      for (int i=0; i<6; i++)
         ccs[i] = i ;

      for (int a=0; a<My_; a++){ // intialize the edge orbits to the solved state
         xorbit xo ;
         uclear(xo) ;
         if (src_ == 'R') // random-state: scramble the wing orbits too
            shuffle(xo.a, xo.a+24, rng) ;
         wingorbits.push_back(xo) ;
      }

      if (src_ == 'U' || alwaysinmem) {
         for (int i=0; i<My; i++) {
            rows[i].resize(Mx) ;
            if (src == 'I' || src == 'C') { // inputmoves, or clear
               for (auto &xo: rows[i])
                  oclear(xo) ;
            } else if (src == 'U') { // orbits have numbers 0..23
               for (auto &xo: rows[i])
                  uclear(xo) ;
            } else if (src == 'R') { // random
               for (auto &xo: rows[i])
                  scramble(xo) ;
            }
         }
      }
      if (src == 'I')
         for (auto mv: moveset)
            domove(mv, *this) ;
   }
   void normalizefacerots() {
      for (int f=0; f<6; f++) {
         if (facerots[f] == 0)
            continue ;
         for (auto &row: rows)
            for (auto &xo: row) {
               switch (facerots[f]) {
                  case 1: xo.rot4(4*f, 4*f+1, 4*f+2, 4*f+3) ; break ;
                  case 2: xo.rot2(4*f, 4*f+1, 4*f+2, 4*f+3) ; break ;
                  case 3: xo.rot4i(4*f, 4*f+1, 4*f+2, 4*f+3) ; break ;
               }
            }
         facerots[f] = 0 ;
      }
   }
   char get(int f, int x, int y) const {
      int ox, oy, oz ;
      coordtospeffz(Mx, My, x, y, ox, oy, oz) ; 
      if (oz == 4)
         return '0'+ccs[f] ;
      if (rows[oy-1].size() == 0)
         return '?' ;
      return (char)('0' + rows[oy-1][ox-1].a[4*f+((oz+facerots[f])&3)]) ;
   }

   void wings_cycle(int p1, int p2, int p3, int p4, int dep){
      auto t = wingorbits[dep].a[p1];
      wingorbits[dep].a[p1] = wingorbits[dep].a[p4];
      wingorbits[dep].a[p4] = wingorbits[dep].a[p3];
      wingorbits[dep].a[p3] = wingorbits[dep].a[p2];
      wingorbits[dep].a[p2] = t;
   }
   void wings_cycle2(int p1, int p2, int p3, int p4, int dep){ 
      swap(wingorbits[dep].a[p1], wingorbits[dep].a[p3]);
      swap(wingorbits[dep].a[p2], wingorbits[dep].a[p4]);
   }

   void rot4hor(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      if (rows[dep].empty()) // centers not materialized in streaming mode; skip
         return ;
      for (int i=0; i<Mx; i++) {
         uchar &c1 = rows[dep][i].a[4*f1+((o1+facerots[f1])&3)] ;
         uchar &c2 = rows[dep][i].a[4*f2+((o2+facerots[f2])&3)] ;
         uchar &c3 = rows[dep][i].a[4*f3+((o3+facerots[f3])&3)] ;
         uchar &c4 = rows[dep][i].a[4*f4+((o4+facerots[f4])&3)] ;
         auto t = c1 ; c1 = c2 ; c2 = c3 ; c3 = c4 ; c4 = t ;
      }
   }
   void rot4ver(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      for (int i=0; i<My; i++) {
         if (rows[i].empty())
            continue ;
         uchar &c1 = rows[i][dep].a[4*f1+((o1+facerots[f1])&3)] ;
         uchar &c2 = rows[i][dep].a[4*f2+((o2+facerots[f2])&3)] ;
         uchar &c3 = rows[i][dep].a[4*f3+((o3+facerots[f3])&3)] ;
         uchar &c4 = rows[i][dep].a[4*f4+((o4+facerots[f4])&3)] ;
         auto t = c1 ; c1 = c2 ; c2 = c3 ; c3 = c4 ; c4 = t ;
      }
   }
   void rot4ver2(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      rot4ver(f1, o1, f2, o2, f3, o3, f4, o4, dep) ;
      rot4ver(f1, o1^2, f2, o2^2, f3, o3^2, f4, o4^2, dep) ;
      uchar &c1 = ccs[f1] ;
      uchar &c2 = ccs[f2] ;
      uchar &c3 = ccs[f3] ;
      uchar &c4 = ccs[f4] ;
      auto t = c1 ; c1 = c2 ; c2 = c3 ; c3 = c4 ; c4 = t ;
   }
   void rot2hor(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      if (rows[dep].empty())
         return ;
      for (int i=0; i<Mx; i++) {
         uchar &c1 = rows[dep][i].a[4*f1+((o1+facerots[f1])&3)] ;
         uchar &c2 = rows[dep][i].a[4*f2+((o2+facerots[f2])&3)] ;
         uchar &c3 = rows[dep][i].a[4*f3+((o3+facerots[f3])&3)] ;
         uchar &c4 = rows[dep][i].a[4*f4+((o4+facerots[f4])&3)] ;
         swap(c1, c3) ;
         swap(c2, c4) ;
      }
   }
   void rot2ver(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      for (int i=0; i<My; i++) {
         if (rows[i].empty())
            continue ;
         uchar &c1 = rows[i][dep].a[4*f1+((o1+facerots[f1])&3)] ;
         uchar &c2 = rows[i][dep].a[4*f2+((o2+facerots[f2])&3)] ;
         uchar &c3 = rows[i][dep].a[4*f3+((o3+facerots[f3])&3)] ;
         uchar &c4 = rows[i][dep].a[4*f4+((o4+facerots[f4])&3)] ;
         swap(c1, c3) ;
         swap(c2, c4) ;
      }
   }
   void rot2ver2(int f1, int o1, int f2, int o2, int f3, int o3, int f4, int o4, int dep) {
      rot2ver(f1, o1, f2, o2, f3, o3, f4, o4, dep) ;
      rot2ver(f1, o1^2, f2, o2^2, f3, o3^2, f4, o4^2, dep) ;
      uchar &c1 = ccs[f1] ;
      uchar &c2 = ccs[f2] ;
      uchar &c3 = ccs[f3] ;
      uchar &c4 = ccs[f4] ;
      swap(c1, c3) ;
      swap(c2, c4) ;
   }
   void getrow(int i) {
      if (alwaysinmem)
         return ;
      if (src != 'R')
         error("! only support random right now for not in memory") ;
      rows[i].resize(Mx) ;
      for (auto &orb: rows[i])
         scramble(orb) ;
   }
   void checkrow(int i) const {
      if ((ll)rows[i].size() != Mx)
         error("! checking a nonexistent row?") ;
      for (auto &orb: rows[i])
         for (int j=0; j<24; j++)
            if (orb.a[j] != (j>>2)) {
               show() ;
               error("! cube is not solved") ;
            }
   }
   void checkrelease(int i) {
      if (alwaysinmem)
         return ;
      checkrow(i) ;
      vector<xorbit> empty ;
      swap(empty, rows[i]) ;
   }

   int wing_color(int oz_wings, int f, int j, int N) {
      static constexpr int colors_of_wings[][2] = { // speffz order on wings
         {0,4}, {0,3}, {0,2}, {0,1}, {1,0}, {1,2}, {1,5}, {1,4},
         {2,0}, {2,3}, {2,5}, {2,1}, {3,0}, {3,4}, {3,5}, {3,2},
         {4,0}, {4,1}, {4,5}, {4,3}, {5,2}, {5,3}, {5,4}, {5,1}
      } ;

      static constexpr int PARTNER[24] = {16,12,8,4,  3,11,23,17, 2,15,20,5,  1,19,21,9, 0,7,22,13, 10,14,18,6};

      // orbit depth k from position j (1..N)
      int k = (j <= N/2) ? (j-1) : (N-j) ;
      int slot = wingorbits[k].a[4*f + oz_wings] ; // 0...23

      if ((j <= N/2)) return colors_of_wings[slot][0] ;
      else return colors_of_wings[ wingorbits[k].a[ PARTNER[4*f + oz_wings] ] ][1]; // PARTNER takes a slot index, not a piece value
   }

   void show_wings() const {
      for (int i=0; i<(int)wingorbits.size(); i++) {
         cout << "edge orbit index " << i+1 << ": [" ;
         for (int j=0; j<24; j++) {
            if (j) cout << "," ;
            cout << (int)wingorbits[i].a[j] ;
         }
         cout << "]" << endl ;
      }
   }
   vector<vector<int>> get_wings(int i){
      vector<vector<int>> r ;
      vector<int> orbit ;
      for (int j=0; j<24; j++)
         orbit.push_back(wingorbits[i].a[j]) ;
      r.push_back(orbit) ;
      return r ;
   }

   cubesizetype Mx, My ;
   vector<vector<xorbit>> rows ;
   vector<xorbit> wingorbits ;
   threecube tc ;
   void show() const ;
   int facerots[6] ;
   uchar ccs[6] ;
   char src ;
} ;
void border() {
   cout << "  " ;
}
void xcube::show() const {
   if (Mx+My > 1000) {
      cout << "Not showing a cube that big" << endl ;
      return ;
   }
   for (int y=1; y<=Mx+My; y++) {
      border() ;
      for (int x=1; x<=Mx+My; x++) {
         cout << " " ;
      }
      border() ;
      for (int x=1; x<=Mx+My; x++)
         cout << get(0, x, y) ;
      cout << endl ;
   }
   cout << endl ;
   for (int y=1; y<=Mx+My; y++) {
      for (int f=1; f<=4; f++) {
         border() ;
         for (int x=1; x<=Mx+My; x++)
            cout << get(f, x, y) ;
      }
      cout << endl ;
   }
   cout << endl ;
   for (int y=1; y<=Mx+My; y++) {
      border() ;
      for (int x=1; x<=Mx+My; x++) {
         cout << " " ;
      }
      border() ;
      for (int x=1; x<=Mx+My; x++)
         cout << get(5, x, y) ;
      cout << endl ;
   }
   cout << endl ;
}
mv parse_move(const char *& s) {
   mv r ;
   r.dep = 1 ;
   while (*s && *s <= ' ')
      s++ ;
   if ('0' <= *s && *s <= '9') {
      r.dep = 0 ;
      while ('0' <= *s && *s <= '9')
         r.dep = 10 * r.dep + *s++ - '0' ;
   }
   auto found = index(faceorder, *s) ;
   if (found == 0) {
      cerr << s << endl ;
      error("! could not parse move") ;
   }
   r.face = found - faceorder ;
   s++ ;
   r.twist = 1 ;
   if (*s && *s > ' ') {
      if ('1' <= *s && *s <= '3')
         r.twist = *s - '0' ;
      else if (*s == '\'')
         r.twist = 3 ;
      else {
         cerr << s << endl ;
         error("! could not parse move (2)") ;
      }
      s++ ;
   }
   while (*s && *s <= ' ')
      s++ ;
   return r ;
}
void readmoveset() {
   string s ;
   while (cin >> s) {
      const char *mstr = s.c_str() ;
      moveset.push_back(parse_move(mstr)) ;
   }
 cout << "Read move set; got " << moveset.size() << endl ;
}
vector<mv> parse_moves(const string &alg) {
   vector<mv> r ;
   const char *s = alg.c_str() ;
   while (*s) {
      while (*s && *s <= ' ')
         s++ ;
      if (*s == 0)
         break ;
      r.push_back(parse_move(s)) ;
   }
   return r ;
}
// we only merge face moves (conjugate moves)
int trymerge(const mv &a, const mv &b, mv &res) {
   if (a.face == b.face && a.dep == b.dep) {
      int tw = (a.twist + b.twist) & 3 ;
      if (tw == 0)
         return 0 ;
      res = a ;
      res.twist = tw ;
      return 1 ;
   } else {
      return 2 ;
   }
}
/*
 *   We retain the sequence of previous face moves so we can try to
 *   cancel them.
 */
vector<mv> prevmvs ;
void trycancellation(mv &m) {
   if (prevmvs.size() == 0) {
      prevmvs.push_back(m) ;
      return ;
   }
   /*
    *   We don't cancel non-face moves at the moment.  We'd have to
    *   be a little smarter in doblock() if we did, because we are
    *   looking at the algorithm template there, not the actual moves
    *   that would be executed.
    */
   if (m.dep > 1) {
      prevmvs.clear() ;
      return ;
   }
   mv merged ;
   int r = trymerge(prevmvs[prevmvs.size()-1], m, merged) ;
   if (r == 2) { // no cancellation
      prevmvs.push_back(m) ;
   } else if (r == 1) {
      cancellations++ ;
      prevmvs[prevmvs.size()-1] = merged ;
   } else {
      cancellations += 2 ;
      prevmvs.pop_back() ;
   }
}
void showmove(const struct mv &mv) {
   if (mv.dep > 1)
      cout << mv.dep ;
   cout << faceorder[mv.face] ;
   if (mv.twist == 2)
      cout << 2 ;
   else if (mv.twist == 3)
      cout << "'" ;
}
void random_move_scramble(cubesizetype N, ll length) {
   moveset.clear() ;
   int lastface = -1 ;
   cout << "SCRAMBLE:" ;
   for (ll i=0; i<length; i++) {
      mv m ;
      do { m.face = rng() % 6 ; } while (m.face == lastface) ; // avoid the identical face twice in a row
      lastface = m.face ;
      m.dep   = 1 + rng() % ((N + 1) / 2) ; 
      m.twist = 1 + rng() % 3 ;
      moveset.push_back(m) ;
      cout << " " ;
      showmove(m) ;
   }
   cout << endl ;
}
// When non-null, every normalized move passed to domove is appended here.
// Used to capture prep moves (corner fix, parity fix, wing cycles) for streaming replay.
vector<mv> *g_prep_recorder = nullptr ;

void domove(mv mv, xcube &cube, bool suppress_cancellations) {
   if (mv.dep > cube.Mx + cube.My + 3 - mv.dep) {
      mv.dep = cube.Mx + cube.My + 3 - mv.dep ;
      mv.face = oppface[mv.face] ;
      mv.twist = (4 - mv.twist) & 3 ;
   }
   if (g_prep_recorder)
      g_prep_recorder->push_back(mv) ;
   assert(mv.dep >= 1) ;
   assert(mv.dep <= cube.Mx + 1) ;
   if (!suppress_cancellations) trycancellation(mv) ;
   if (mv.dep == 1) { // we can't easily ignore wings on face moves
      cube.tc.domove(mv) ;
      cube.facerots[mv.face] = (cube.facerots[mv.face] + 4 - mv.twist) & 3 ;
      for (int m=0; m<cube.My; m++){
         for (int i=0; i<mv.twist; i++) {
            if (mv.face == 0) { // U
               cube.wings_cycle(0,1,2,3, m);
               cube.wings_cycle(16,12,8,4, m);
            } else if (mv.face == 1) { // L
               cube.wings_cycle(4,5,6,7, m);
               cube.wings_cycle(3,11,23,17, m);
            } else if (mv.face == 2) { // F
               cube.wings_cycle(8,9,10,11, m);
               cube.wings_cycle(2,15,20,5, m);
            } else if (mv.face == 3) { // R
               cube.wings_cycle(12,13,14,15, m);
               cube.wings_cycle(1,19,21,9, m);
            } else if (mv.face == 4) { // B
               cube.wings_cycle(16,17,18,19, m);
               cube.wings_cycle(0,7,22,13, m);
            } else if (mv.face == 5) { // D
               cube.wings_cycle(20,21,22,23, m);
               cube.wings_cycle(10,14,18,6, m);
            }
         }
      }
      return ;
   }
   assert(mv.twist > 0) ;
   assert(mv.twist < 4) ;
   int dep = mv.dep - 2 ;
   /*
    *   For odd cubes, when doing the exact middle slice, things are
    *   different; this is ugly but okay.
    */
   if (dep >= cube.My) {
      cube.tc.domove({2, mv.face, mv.twist}) ;
      switch (mv.face * 3 + mv.twist) {
case 1: // U
         cube.rot4ver2(1, 1, 2, 1, 3, 1, 4, 1, dep) ;
         break ;
case 2: // U2
         cube.rot2ver2(1, 1, 2, 1, 3, 1, 4, 1, dep) ;
         break ;
case 3: // U'
         cube.rot4ver2(4, 1, 3, 1, 2, 1, 1, 1, dep) ;
         break ;
case 4: // L
         cube.rot4ver2(2, 0, 0, 0, 4, 2, 5, 0, dep) ;
         break ;
case 5: // L2
         cube.rot2ver2(2, 0, 5, 0, 4, 2, 0, 0, dep) ;
         break ;
case 6: // L'
         cube.rot4ver2(2, 0, 5, 0, 4, 2, 0, 0, dep) ;
         break ;
case 7: // F
         cube.rot4ver2(0, 3, 1, 2, 5, 1, 3, 0, dep) ;
         break ;
case 8: // F2
         cube.rot2ver2(0, 3, 3, 0, 5, 1, 1, 2, dep) ;
         break ;
case 9: // F'
         cube.rot4ver2(0, 3, 3, 0, 5, 1, 1, 2, dep) ;
         break ;
case 10: // R
         cube.rot4ver2(2, 2, 5, 2, 4, 0, 0, 2, dep) ;
         break ;
case 11: // R2
         cube.rot2ver2(2, 2, 5, 2, 4, 0, 0, 2, dep) ;
         break ;
case 12: // R'
         cube.rot4ver2(2, 2, 0, 2, 4, 0, 5, 2, dep) ;
         break ;
case 13: // B
         cube.rot4ver2(0, 1, 3, 2, 5, 3, 1, 0, dep) ;
         break ;
case 14: // B2
         cube.rot2ver2(0, 1, 3, 2, 5, 3, 1, 0, dep) ;
         break ;
case 15: // B'
         cube.rot4ver2(0, 1, 1, 0, 5, 3, 3, 2, dep) ;
         break ;
case 16: // D
         cube.rot4ver2(4, 3, 3, 3, 2, 3, 1, 3, dep) ;
         break ;
case 17: // D2
         cube.rot2ver2(1, 3, 2, 3, 3, 3, 4, 3, dep) ;
         break ;
case 18: // D'
         cube.rot4ver2(4, 3, 1, 3, 2, 3, 3, 3, dep) ;
         break ;
      }
   } else {
      switch (mv.face * 3 + mv.twist) {
case 1: // U
         cube.rot4hor(1, 0, 2, 0, 3, 0, 4, 0, dep) ;
         cube.rot4ver(1, 1, 2, 1, 3, 1, 4, 1, dep) ;
         cube.wings_cycle(5,17,13,9, dep);
         break ;
case 2: // U2
         cube.rot2hor(1, 0, 2, 0, 3, 0, 4, 0, dep) ;
         cube.rot2ver(1, 1, 2, 1, 3, 1, 4, 1, dep) ;
         cube.wings_cycle2(5,17,13,9, dep);
         break ;
case 3: // U'
         cube.rot4hor(4, 0, 3, 0, 2, 0, 1, 0, dep) ;
         cube.rot4ver(4, 1, 3, 1, 2, 1, 1, 1, dep) ;
         cube.wings_cycle(5,9,13,17, dep);
         break ;
case 4: // L
         cube.rot4hor(2, 3, 0, 3, 4, 1, 5, 3, dep) ;
         cube.rot4ver(2, 0, 0, 0, 4, 2, 5, 0, dep) ;
         cube.wings_cycle(0,8,20,18, dep);
         break ;
case 5: // L2
         cube.rot2hor(2, 3, 5, 3, 4, 1, 0, 3, dep) ;
         cube.rot2ver(2, 0, 5, 0, 4, 2, 0, 0, dep) ;
         cube.wings_cycle2(0,8,20,18, dep);
         break ;
case 6: // L'
         cube.rot4hor(2, 3, 5, 3, 4, 1, 0, 3, dep) ;
         cube.rot4ver(2, 0, 5, 0, 4, 2, 0, 0, dep) ;
         cube.wings_cycle(0,18,20,8, dep);
         break ;
case 7: // F
         cube.rot4hor(0, 2, 1, 1, 5, 0, 3, 3, dep) ;
         cube.rot4ver(0, 3, 1, 2, 5, 1, 3, 0, dep) ;
         cube.wings_cycle(3,12,21,6, dep);
         break ;
case 8: // F2
         cube.rot2hor(0, 2, 3, 3, 5, 0, 1, 1, dep) ;
         cube.rot2ver(0, 3, 3, 0, 5, 1, 1, 2, dep) ;
         cube.wings_cycle2(3,12,21,6, dep);
         break ;
case 9: // F'
         cube.rot4hor(0, 2, 3, 3, 5, 0, 1, 1, dep) ;
         cube.rot4ver(0, 3, 3, 0, 5, 1, 1, 2, dep) ;
         cube.wings_cycle(3,6,21,12, dep);
         break ;
case 10: // R
         cube.rot4hor(2, 1, 5, 1, 4, 3, 0, 1, dep) ;
         cube.rot4ver(2, 2, 5, 2, 4, 0, 0, 2, dep) ;
         cube.wings_cycle(2, 16, 22, 10, dep);
         break ;
case 11: // R2
         cube.rot2hor(2, 1, 5, 1, 4, 3, 0, 1, dep) ;
         cube.rot2ver(2, 2, 5, 2, 4, 0, 0, 2, dep) ;
         cube.wings_cycle2(2, 16, 22, 10, dep);
         break ;
case 12: // R'
         cube.rot4hor(2, 1, 0, 1, 4, 3, 5, 1, dep) ;
         cube.rot4ver(2, 2, 0, 2, 4, 0, 5, 2, dep) ;
         cube.wings_cycle(2, 10, 22, 16, dep);  
         break ;
case 13: // B
         cube.rot4hor(0, 0, 3, 1, 5, 2, 1, 3, dep) ;
         cube.rot4ver(0, 1, 3, 2, 5, 3, 1, 0, dep) ;
         cube.wings_cycle(1, 4, 23, 14, dep);
         break ;
case 14: // B2
         cube.rot2hor(0, 0, 3, 1, 5, 2, 1, 3, dep) ;
         cube.rot2ver(0, 1, 3, 2, 5, 3, 1, 0, dep) ;
         cube.wings_cycle2(1, 4, 23, 14, dep);
         break ;
case 15: // B'
         cube.rot4hor(0, 0, 1, 3, 5, 2, 3, 1, dep) ;
         cube.rot4ver(0, 1, 1, 0, 5, 3, 3, 2, dep) ;
         cube.wings_cycle(1, 14, 23, 4, dep);
         break ;
case 16: // D
         cube.rot4hor(4, 2, 3, 2, 2, 2, 1, 2, dep) ;
         cube.rot4ver(4, 3, 3, 3, 2, 3, 1, 3, dep) ;
         cube.wings_cycle(7, 11, 15, 19, dep);   //### was (7,19,15,11): wrong direction
         break ;
case 17: // D2
         cube.rot2hor(1, 2, 2, 2, 3, 2, 4, 2, dep) ;
         cube.rot2ver(1, 3, 2, 3, 3, 3, 4, 3, dep) ;
         cube.wings_cycle2(7, 19, 15, 11, dep);
         break ;
case 18: // D'
         cube.rot4hor(4, 2, 1, 2, 2, 2, 3, 2, dep) ;
         cube.rot4ver(4, 3, 1, 3, 2, 3, 3, 3, dep) ;
         cube.wings_cycle(7, 19, 15, 11, dep);   //### was (7,11,15,19): wrong direction
         break ;
      }
   }
}
#if WRITEFRAMES
int frameseq = 0 ;
uchar *framebuf = 0 ;
ll framew, frameh ;
int facexoff[] = { 1, 0, 1, 2, 3, 1 } ;
int faceyoff[] = { 0, 1, 1, 1, 1, 2 } ;
/*               w    o     g    r    b    y  grey dark_orange*/
int colorr[] = { 255, 255,   0, 255,   0, 255,  128, 255} ;
int colorg[] = { 255, 165, 160,   0,   0, 255,  128, 128} ;
int colorb[] = { 255,   0,   0,   0, 255,   0,  128, 0} ;
const int S = 24 ;
void placepixel(int col, int x, int y) {
   int st = (S >= 5 ? 1 : 0) ;
   int r=colorr[col], g=colorg[col], b=colorb[col] ;
   for (int dy=st; dy<S; dy++) for (int dx=st; dx<S; dx++) {
      ll o = ((ll)(y*S+dy)*framew + (x*S+dx)) * 3 ;
      framebuf[o]=r; framebuf[o+1]=g; framebuf[o+2]=b;
   }
}
void writeframe(xcube &cube) {
   writemovetarget += writemoveinc ;
   if (cube.Mx + cube.My > 1000000000)
      error("! cube too large for frame output") ;
   int N = (int)(cube.Mx + cube.My) ;
   if (framebuf == 0) {
      framew = 4LL * (N+2) * S ;
      frameh = 3LL * (N+2) * S ;
      ll framemem = 3 * framew * frameh ;
      framebuf = (uchar *)calloc(1, framemem) ;
      if (framebuf == 0)
         error("! no memory for the frame buffer") ;
   }
   for (int i=0; i<cube.My; i++) {
      auto &row = cube.rows[i] ;
      for (int j=0; j<cube.Mx; j++) {
         auto &xo = row[j] ;
         for (int f=0; f<6; f++) {
            int fx = facexoff[f] * (N+2) +1 ;
            int fy = faceyoff[f] * (N+2) +1 ;
            placepixel(xo.a[4*f+((cube.facerots[f])&3)], fx+j, fy+i) ;
            placepixel(xo.a[4*f+((1+cube.facerots[f])&3)], fx+N-1-i, fy+j) ;
            placepixel(xo.a[4*f+((2+cube.facerots[f])&3)], fx+N-1-j, fy+N-1-i) ;
            placepixel(xo.a[4*f+((3+cube.facerots[f])&3)], fx+i, fy+N-1-j) ;
         }
      }
   }
   for (int f=0; f<6; f++) //###
      for (int i=1; i<=N; i++) { // skip corners
         if (i+i==N+1) // skip middle edges
            continue ;
         placepixel(cube.wing_color(0, f, i, N), facexoff[f]*(N+2)+i,        faceyoff[f]*(N+2)) ;
         placepixel(cube.wing_color(1, f, i, N), facexoff[f]*(N+2)+N+1,      faceyoff[f]*(N+2)+i) ;
         placepixel(cube.wing_color(2, f, i, N), facexoff[f]*(N+2)+(N+1-i),  faceyoff[f]*(N+2)+N+1) ;
         placepixel(cube.wing_color(3, f, i, N), facexoff[f]*(N+2),          faceyoff[f]*(N+2)+(N+1-i)) ;
      }
   auto &threecube = cube.tc ;
   for (int f=0; f<6; f++) {
      int fx = facexoff[f] * (N+2);
      int fy = faceyoff[f] * (N+2);
      placepixel(threecube.corner.a[4*f], fx, fy) ;
      placepixel(threecube.corner.a[4*f+1], fx+N+1, fy) ;
      placepixel(threecube.corner.a[4*f+2], fx+N+1, fy+N+1) ;
      placepixel(threecube.corner.a[4*f+3], fx, fy+N+1) ;
      if (N & 1) {
         int half = (N+1)/2 ;
         placepixel(threecube.edge.a[4*f], fx+half, fy) ;
         placepixel(threecube.edge.a[4*f+1], fx+N+1, fy+half) ;
         placepixel(threecube.edge.a[4*f+2], fx+half, fy+N+1) ;
         placepixel(threecube.edge.a[4*f+3], fx, fy+half) ;
         placepixel(threecube.centercenters[f], fx+half, fy+half) ;
      }
   }
   string framename = "frame-" ;
   cout << " " << frameseq << flush ;
   for (int dig=1000000; dig>=1; dig /= 10)
      framename.push_back((char)(frameseq / dig % 10 + '0')) ;
   framename += ".png" ;
   frameseq++ ;
   int t = fpng::fpng_encode_image_to_file(framename.c_str(), framebuf,
      framew, frameh, 3) ;
   if (t > 1) {
      cerr << t << endl ;
      error("! problem writing frame") ;
   }
}
#endif
/*
 *   Using elements 1 and 2 as generators, expand the set of rotations.
 *   We can encode each permutation of faces as a 32-bit int, and use
 *   the explicit static array as a bfs queue, so this is pretty easy.
 */
int encoderot(int *p) {
   return p[0] * 6 + p[1] ;
}
void expandrotations() {
   const int MAXENC = 6*6 ;
   char seen[MAXENC] ;
   if (rots[23][0] != 0 || rots[23][1] != 0)
      return ;
   for (int i=0; i<MAXENC; i++)
      seen[i] = 0 ;
   for (int i=0; i<3; i++)
      seen[encoderot(rots[i])] = 1 ;
   int qp = 3 ;
   int w[6] ;
   for (int qg=1; qg<24; qg++) {
      for (int g=1; g<3; g++) {
         for (int i=0; i<6; i++)
            w[i] = rots[qg][rots[g][i]] ;
         int ev = encoderot(w) ;
         if (seen[ev] == 0) {
            seen[ev] = 1 ;
            for (int i=0; i<6; i++)
               rots[qp][i] = w[i] ;
            qp++ ;
            if (qp > 24)
               error("! blew up expanding rotations") ;
         }
      }
   }
   if (qp != 24)
      error("! did not find the correct count of rotations") ;
}
/*
 *   How much does this algorithm cost?
 */
int oppositemoves(const mv &a, const mv &b) {
   return (a.dep == b.dep && a.face == b.face && a.twist + b.twist == 4) ;
}
struct blockalgo {
   vector<mv> mvs ;
   xorbit perm ;
   ll cost ; // vhoc each has three base-10 digits
   int v, h, o, c ;
   blockalgo(const vector<mv> &mvs_, const xorbit &perm_) :
                   mvs(mvs_), perm(perm_), cost(0), v(0), h(0), o(0), c(0) {

      int end = mvs.size()-1 ;
      int st = 0 ;
      while(st < end && mvs[st].dep == 1 && oppositemoves(mvs[st], mvs[end])) {
         end-- ;
         st++ ;
         c += 2 ;
      }
      for (int i=st; i<=end; i++) {
         switch (mvs[i].dep) {
case 1:     o++ ; break ;
case 2:     h++ ; break ;
case 3:     v++ ; break ;
default:    error("! unexpected depth in algorithm") ;
         }
      }
      cost = c + 10000 * (o + 10000 * (h + 10000LL * v)) ;
   }
   void show() const {
      for (auto m: mvs) {
         cout << " " ;
         showmove(m) ;
      }
   }
} ;
vector<blockalgo> loadalgs(const vector<string> &algs, int m=0) {
   expandrotations() ;
   xcube cube(4, 5, 'U') ; // get the special version with 0..23
   vector<mv> mappedmoves ;
   vector<blockalgo> r ;
   for (int i=0; i<(int)algs.size(); i++) {
      auto mvs = parse_moves(algs[i]) ;
      uclear(cube.rows[0][1]) ;
      mappedmoves.clear() ;
      for (auto mv: mvs) {
         mv.face = rots[m][mv.face] ;
         mappedmoves.push_back(mv) ;
         domove(mv, cube, true) ; // suppress_cancellations: loadalgs is called from parallel threads
      }
      // Test that algorithms don't affect other orbits
      for (ll j=0; j<cube.My; j++)
         for (ll k=0; k<cube.Mx; k++)
            if (j != 0 && k != 1) { 
               for (int p=0; p<24; p++)
                  if (cube.rows[j][k].a[p] != p) {
                     cout << "Failed alg " << i << " at " << j << " " << k << " touches " << p << endl ;
                     error("! algorithm hurts other orbits?") ;
                  }
            }
      int changed = 0 ;
      for (int i=0; i<24; i++)
         if (cube.rows[0][1].a[i] != i)
            changed++ ;
      if (changed == 0)
         error("! algorithm doesn't change any cubies") ;
      r.push_back(blockalgo(mappedmoves, cube.rows[0][1])) ;
   }
   return r ;
}
/*
 *   We always buffer one block, so we can cancel moves between blocks.
 */
vector<mv> curblock, prevblock ;
void sendblock() {
   while (curblock.size() > 0 && prevblock.size() > 0) {
      mv merged ;
      int r = trymerge(prevblock[prevblock.size()-1], curblock[0], merged) ;
      if (r == 0) { // cancelled completely
         curblock.erase(curblock.begin()) ;
         prevblock.pop_back() ;
      } else if (r == 1) { // merged to just 1
         curblock[0] = merged ;
         prevblock.pop_back() ;
      } else { // nothing doing
         break ;
      }
   }
   if (prevblock.size()) {
      for (auto mv: prevblock) {
         cout << " " ;
         showmove(mv) ;
      }
      cout << endl ;
   }
   swap(prevblock, curblock) ;
   curblock.clear() ;
}
ll doblock(xcube &cube, blockalgo &alg, vector<int> &ys, vector<int> &xs) {
#if WRITEFRAMES
   if (writemoveinc && runningmovecount >= writemovetarget)
      writeframe(cube) ;
#endif
   ll oldcancellations = cancellations ;
   blocks++ ;
   if (writemoves) {
      for (auto mv: alg.mvs) {
         if (mv.dep == 2) {
            for (auto v:ys) {
               mv.dep = v + 2 ;
               curblock.push_back(mv) ;
            }
         } else if (mv.dep == 3) {
            for (auto v:xs) {
               mv.dep = v + 2 ;
               curblock.push_back(mv) ;
            }
         } else {
            curblock.push_back(mv) ;
         }
      }
      sendblock() ;
   }
   if (!fast) {
      for (auto mv: alg.mvs) {
         if (mv.dep == 2) {
            for (auto v:ys) {
               mv.dep = v + 2 ;
               domove(mv, cube) ;
            }
         } else if (mv.dep == 3) {
            for (auto v:xs) {
               mv.dep = v + 2 ;
               domove(mv, cube) ;
            }
         } else {
            domove(mv, cube)  ;
         }
      }
   } else {
      for (int i=0; i<alg.c; i++)
         trycancellation(alg.mvs[i]) ;
      prevmvs.clear() ;
      for (int i=alg.c-1; i>=0; i--)
         trycancellation(alg.mvs[alg.mvs.size()-1-i]) ;
   }
   ll mvs = alg.o + alg.v * xs.size() + alg.h * ys.size() + alg.c -
                                         (cancellations - oldcancellations) ;
   runningmovecount += mvs ;
   return mvs ;
}
struct Stats {
   ll hits, moves ;
   double expectedcost ;
   blockalgo *ba ;
   int w1, w2 ;
   cubesizetype count1, count2 ;
   void showstats(int id, double denom) {
      if (hits) {
         cout << "Phase " << id << " p " << hits/denom << " mvs " << moves << " expected " << expectedcost ;
         if (count1 > 0 || count2 > 0) {
            cout << " widths:" ;
            if (count1 > 0) cout << " " << w1 << "x" << count1 ;
            if (count2 > 0) cout << " " << w2 << "x" << count2 ;
         }
         cout << endl ;
         cout << "   " ; ba->show() ; cout << endl ;
      }
   }
} ;
vector<Stats> byphase ;
const double INFCOST = 1e308 ;
/*
 *   Here we optimize for a discrete w.  We could do it in O(N log N)
 *   by considering each row, but we don't need to; the cost function
 *   is convex in w, so an optimal solution will always have exactly
 *   one or two distinct values, and if there are two, they will differ
 *   by one.  Further, one will dominate.  So we just try all possibilities
 *   of a major w and a minor w.
 *
 *   We also keep the argument information so we can recalculate a plan
 *   for a different My.  We need to do this for the last partition if we
 *   partition unevenly; for instance, if we solve for N=1024 and we
 *   choose p=100, there will be a last part with only 11 rows.
 */
struct wplan {
   double cost ;
   cubesizetype My, count1, count2 ;
   int w1, w2 ;
   vector<double> byw ;
   void plan(cubesizetype My_) {
      if (My == My_)
         return ;
      My = My_ ;
      cost = INFCOST ;
      for (int i=1; i<(int)byw.size(); i++) {
         // use mostly i's but some i-1's to make up the difference
         if (i > 1) {
            cubesizetype mless = (i - My % i) % i ; // this many lessor ones
            cubesizetype mi = (My - mless * (i - 1)) / i ;
            if (mi * i + mless * (i - 1) != My)
               error("! bad calculation") ;
            if (mi >= 0) {
               double t = mless * byw[i-1] + mi * byw[i] ;
               if (t < cost) {
                  cost = t ;
                  w1 = i-1 ;
                  w2 = i ;
                  count1 = mless ;
                  count2 = mi ;
               }
            }
         }
         // use mostly i's but some i+1's to make up the difference
         cubesizetype mgreat = My % i ;
         cubesizetype c1 = My / i - mgreat ;
         if (c1 >= 0 && (mgreat == 0 || i + 1 < (int)byw.size())) {
            double t = c1 * byw[i] ;
            if (mgreat > 0)
               t += mgreat * byw[i+1] ;
            if (t < cost) {
               cost = t ;
               w1 = i ;
               w2 = i+1 ;
               count1 = c1 ;
               count2 = mgreat ;
            }
         }
      }
 if (count1 * w1 + count2 * w2 != My) cout << "count1 " << count1 << " w1 " << w1 << " count2 " << count2 << " w2 " << w2 << " My " << My << endl ;
      assert(count1 * w1 + count2 * w2 == My) ;
   }
} ;
struct netgen ;
struct reductionwplan {
   double cost ;
   netgen *net ;
   vector<wplan> plan ;
} ;
wplan optw(const vector<double> &byw, cubesizetype My) {
   wplan r = {INFCOST, -1, -1LL, -1LL, -1, -1, byw} ;
   r.plan(My) ;
   return r ;
}
wplan avcost(ll cost, cubesizetype Mx, cubesizetype My, double prob) {
   vector<double> byw = {0} ;
   int v = cost / 1000000000000LL ;
   int h = cost / 100000000 % 10000 ;
   int o = cost / 10000 % 10000 ;
   int c = cost % 10000 ;
   double q = 1 - prob ;
   double ratio = prob / q ;
   double qpow = q ;   // tracks q^w
   for (int w=1; ; w++) {
      double t = 0 ;
      double p1 = qpow * ratio ;   // prob^1 * q^(w-1) for bs=1
      double choose = w ;
      for (int bs=1; bs<=w; bs++) {
 //      double any = 1-pow(1-p1, Mx) ;  // note:  might be very inaccurate
         double any = -expm1(Mx * log1p(-p1));  // much better accuracy
         t += choose * (v * p1 * Mx + h * any * bs + o * any) ;
         choose = choose * (w - bs) / (bs + 1) ;
         p1 *= ratio ;
      }
      byw.push_back(t) ;
      if ((Mx >> (w + 1)) == 0)
         break ;
      qpow *= q ;
   }
   auto r = optw(byw, My) ;
   r.cost += c ;
   return r ;
}
wplan wccost(ll cost, cubesizetype Mx, cubesizetype My) {
   vector<double> byw = {0} ;
   int v = cost / 1000000000000LL ;
   int h = cost / 100000000 % 10000 ;
   int o = cost / 10000 % 10000 ;
   int c = cost % 10000 ;
   for (int w=1; ; w++) {
      double t = (v * Mx + h * w * (1LL << w) * 0.5 + o * ((1LL << w) - 1)) ;
      byw.push_back(t) ;
      if ((Mx >> (w + 1)) == 0)
         break ;
   }
   auto r = optw(byw, My) ;
   r.cost += c ;
   return r ;
}
/*
 *   Often we need to canonicalize or remap numbers so they are in order.
 *   This little struct does the work.
 */
struct remapper {
   remapper() {
      for (int i=0; i<12; i++)
         newval[i] = 0 ;
      at = '0' ;
   }
   char get(char c) {
      if (c == 'x')
         return 'x' ;
      if (newval[c-'0'] == 0)
         newval[c-'0'] = at++ ;
      return newval[c-'0'] ;
   }
   int at ;
   char newval[12] ;
   static void remap(string &k) {
      remapper remap ;
      for (int i=0; i<6; i++)
         k[i] = remap.get(k[i]) ;
   }
} ;
/*
 *   This is our general net class for a particular reduction (such as
 *   24/12/12 or 8/4/4), but with specific faces.  Input is a six-character
 *   description of a reduction which is a combination of x characters
 *   and digits 0, 1, and possibly 2.  The characters correspond directly
 *   to the faces in the Speffz face ordering ULFRBD.  An x character means
 *   that face is unchanged.  The digit faces are the ones that we sort,
 *   and 0's, 1's, and optional 2's indicate the groups we sort them into.
 *   So 012xxx means that the U, L, and F faces already have the centers
 *   sorted into just U, L, and F pieces (but possibly scrambled), and
 *   we are sorting them so all U pieces are on the U face, and so on for
 *   the other faces.  A 000111 description means no sorting has been
 *   done yet, but we are placing the U, L, and F pieces onto the U, L, and
 *   F faces (which leaves the R, B, and D pieces on the R, B, and D faces).
 *   So one particular sorting strategy might be 000111 followed by
 *   012xxx and xxx012.  Another might be 010100 followed by x0x1xx and
 *   0x0x11 and 0x1xxx and xxxx01.
 *
 *   We use our rotational intelligence above to derive netgens that are
 *   rotationally equivalent to and from each other.
 *
 *   There may be different netgens with different cost structures; we
 *   have a strategy later that picks the one most appropriate to a
 *   particular cube size.
 */
/*
 *   Pack/unpack tables shared across all netgens with the same digcounts.
 *   Building these is moderately expensive, so we cache one copy per unique
 *   digcount tuple and hand out shared_ptr<const PackTables> to each netgen.
 */
struct PackTables {
   vector<int> packhi, packlo ;
   vector<ll> unpack ;
   vector<ll> basepow ;
   int halfmod, goal, diggoal, base ;
   PackTables(int digcounts[4]) {
      base = (digcounts[3] ? 4 : digcounts[2] ? 3 : 2) ;
      int nonx = digcounts[0] + digcounts[1] + digcounts[2] + digcounts[3] ;
      basepow.resize(4 * nonx) ;
      basepow[0] = 1 ;
      for (int i=1; i<4*nonx; i++)
         basepow[i] = basepow[i-1] * base ;
      ll cases = 1, numer = 4 * nonx ;
      for (int i=0; i<4; i++)
         for (int j=0; j<4*digcounts[i]; j++) {
            cases = cases * numer / (j + 1) ;
            if (cases > 1000000000000000LL)
               error("! bad key; too many cases") ;
            numer-- ;
         }
      if (cases > 100000000)
         error("! bad key; too many cases") ;
      goal = cases ;
      halfmod = 1 ;
      for (int i=0; i<4*nonx; i += 2)
         halfmod *= base ;
      int digval[] = {1, 100, 10000, 1000000} ;
      vector<int> digsums(halfmod) ;
      digsums[0] = 4 * nonx / 2 ;
      for (int i=1; i<halfmod; i++)
         digsums[i] = digval[i % base] + digsums[i / base] - 1 ;
      diggoal = 4 * (digval[0] * digcounts[0] + digval[1] * digcounts[1] +
                     digval[2] * digcounts[2] + digval[3] * digcounts[3]) ;
      int at = 0 ;
      packhi.resize(halfmod, -1) ;
      packlo.resize(halfmod, -1) ;
      vector<pair<int,int>> sortme(halfmod) ;
      for (int th=0; th<halfmod; th++)
         sortme[th] = {digsums[th], th} ;
      sort(sortme.begin(), sortme.end()) ;
      int bw = halfmod - 1 ;
      for (int thi=0; thi<halfmod; thi++) {
         while (bw >= 0 && sortme[bw].first + sortme[thi].first >= diggoal)
            bw-- ;
         for (int bhi=bw+1; bhi<halfmod && sortme[bhi].first + sortme[thi].first == diggoal; bhi++) {
            if (packhi[sortme[thi].second] == -1)
               packhi[sortme[thi].second] = at ;
            if (packlo[sortme[bhi].second] == -1)
               packlo[sortme[bhi].second] = at - packhi[sortme[thi].second] ;
            unpack.push_back(sortme[thi].second * (ll)halfmod + sortme[bhi].second) ;
            at++ ;
         }
      }
      if (goal != at) {
         cout << "Digvals " << digval[0] << " " << digval[1] << " " << digval[2] << endl ;
         cout << "Digcounts " << digcounts[0] << " " << digcounts[1] << " " << digcounts[2] << endl ;
         cout << "Goal " << goal << " at " << at << " diggoal " << diggoal << endl ;
         error("! failed in goal computation") ;
      }
   }
} ;
static mutex pack_cache_mu ;
static map<array<int,4>, shared_ptr<const PackTables>> pack_cache ;
static shared_ptr<const PackTables> get_pack_tables(int digcounts[4]) {
   array<int,4> key = {digcounts[0], digcounts[1], digcounts[2], digcounts[3]} ;
   lock_guard<mutex> lk(pack_cache_mu) ;
   auto it = pack_cache.find(key) ;
   if (it != pack_cache.end())
      return it->second ;
   auto pt = make_shared<PackTables>(digcounts) ;
   return pack_cache.emplace(key, std::move(pt)).first->second ;
}
struct netgen {
   const static int MINWRITE = 50 ;
   int STAGE ;
   static const int n ;
   shared_ptr<const PackTables> pt ;
   vector<double> stageprob ;
   string key, filename ;
   vector<blockalgo> algs ;
   vector<ll> algcosts ;
   vector<string> stringalgs ;
   int halfmod, goal, diggoal ;
   int base ; // 2, 3, or 4
   int tolowdig[24] ;
   ll enccost ;
   int m ;
   bool optimal ;
   bool bad = false ;
   int used_perms = -1 ;
   int max_activations_seen = -1 ;
   netgen(string key_, const vector<string> &stringalgs_, const string &filename_, int m_=0, bool optimal_=false) {
      optimal = optimal_ ;
      /*
       *   Pick apart the key, build tolowdig, and fetch shared pack tables.
       */
      m = m_ ;
      key.resize(6, 0) ;
      for (int i=0; i<6; i++)
         key[rots[m][i]] = key_[i] ;
      filename = filename_ ;
      stringalgs = stringalgs_ ;
      STAGE = stringalgs.size() ;
      int digcounts[4] = {} ;
      if (key.size() != 6)
         error("! invalid key 1") ;
      int at = 0 ;
      int f = 0 ;
      for (auto c: key) {
         if (c == 'x') {
            for (int i=0; i<4; i++)
               tolowdig[4*f+i] = -1 ;
         } else if ('0' <= c && c <= '3') {
            digcounts[c-'0']++ ;
            for (int i=0; i<4; i++)
               tolowdig[4*f+i] = at++ ;
         } else {
            cerr << key << endl ;
            error("! bad char in netgen key") ;
         }
         f++ ;
      }
      if (digcounts[0] == 0 || digcounts[1] == 0)
         error("! invalid key 2") ;
      pt = get_pack_tables(digcounts) ;
      halfmod = pt->halfmod ;
      goal    = pt->goal ;
      diggoal = pt->diggoal ;
      base    = pt->base ;
      enccost = 0 ;
   }
   void ensure_algs_loaded() {
      if (!algs.empty()) return ;
      algs = loadalgs(stringalgs, m) ;
      enccost = 0 ;
      for (auto &alg: algs) {
         enccost += alg.cost ;
         for (int f=0; f<6; f++)
            if (key[f] == 'x')
               for (int i=0; i<4; i++)
                  if (alg.perm.a[4*f+i] != 4*f+i)
                     error("! an algorithm in netgen modifies the wrong face") ;
      }
      if (algcosts.empty())
         for (auto &alg: algs)
            algcosts.push_back(alg.cost) ;
   }
   /*
    *   This is kind of slow.  We can speed it up if we want, later.
    */
   ll permdigs(ll v, xorbit &xo) {
      ll r = 0 ;
      if (base == 2 && halfmod == 4096) { // assume this gets unrolled
         for (int i=0; i<24; i++) {
            r += (v & 1) << xo.a[i] ;
            v >>= 1 ;
         }
      } else {
         for (int i=0; i<24; i += 4) {
            if (tolowdig[i] >= 0) {
               r += pt->basepow[tolowdig[xo.a[i]]] * (v % base) ;
               v /= base ;
               r += pt->basepow[tolowdig[xo.a[i+1]]] * (v % base) ;
               v /= base ;
               r += pt->basepow[tolowdig[xo.a[i+2]]] * (v % base) ;
               v /= base ;
               r += pt->basepow[tolowdig[xo.a[i+3]]] * (v % base) ;
               v /= base ;
            }
         }
      }
      return r ;
   }
/*
 *   Precompute per-face lookup tables for fast plan passes.
 *   tab[g][val] is the k-contribution from the g-th non-x face when its
 *   4 stickers encode to 'val' in little-endian base^4 (= base*base*base*base).
 *   Returns the number of non-x face groups (= number of non-x faces).
 */
   int build_fast_tab(const xorbit &xo, ll tab[][256]) const {
      const int base4 = base * base * base * base ;
      int ng = 0 ;
      for (int f = 0; f < 6; f++) {
         if (tolowdig[4*f] < 0) continue ;
         for (int val = 0; val < base4; val++) {
            ll contrib = 0 ; int tmp = val ;
            for (int di = 0; di < 4; di++) {
               contrib += pt->basepow[tolowdig[xo.a[4*f + di]]] * (tmp % base) ;
               tmp /= base ;
            }
            tab[ng][val] = contrib ;
         }
         ng++ ;
      }
      return ng ;
   }
/*
 *   Plan how to execute this based on the permutations from the
 *   algorithms.
 */
   vector<planmask> plans ;
   int starti_for_plan() const {
      ll start = 0 ;
      for (int i=5; i>=0; i--)
         if (key[i] != 'x')
            for (int j=0; j<4; j++)
               start = start * base + key[i] - '0' ;
      return pt->packhi[start/halfmod] + pt->packlo[start%halfmod] ;
   }
   void finish_plan() {
      for (int i=0; i<goal; i++)
         if (plans[i].is_none()) {
            cerr << "Inconsistency in file: " << filename << endl ;
            cerr << "  key=" << key << " goal=" << goal << " STAGE=" << STAGE << " optimal=" << optimal << endl ;
            cerr << "  Could not get to state " << i << endl ;
            bad = true ;
            return ;
         }
      vector<int> cnts(STAGE) ;
      for (int i=0; i<goal; i++)
         for (int j=0; j<STAGE; j++)
            if (plans[i].test(j))
               cnts[j]++ ;
      used_perms = 0 ;
      for (int i=0; i<STAGE; i++)
         if (cnts[i] > 0)
            used_perms++ ;
      max_activations_seen = 0 ;
      for (int i=0; i<goal; i++) {
         int pc = plans[i].popcount() ;
         if (pc > max_activations_seen)
            max_activations_seen = pc ;
      }
      stageprob.clear() ;
      for (int i=0; i<STAGE; i++)
         stageprob.push_back(cnts[i]/(double)goal) ;
   }
   void plan_greedy() {
      ensure_algs_loaded() ;
      double t0 = walltime() ;
      if (pt->unpack.empty())
         error("! unpack should not be empty") ;
      vector<planmask> next(goal, planmask::none()) ;
      vector<ll> nextcost(goal, LLONG_MAX) ;
      vector<ll> curcost(goal, LLONG_MAX) ;
      plans.resize(goal, planmask::none()) ;
      int starti = starti_for_plan() ;
      plans[starti] = planmask{} ;
      curcost[starti] = 0 ;
      int covered = 1 ;
      const bool pow2base = !(base & (base - 1)) ;
      const int hmshift = pow2base ? __builtin_ctz(halfmod) : 0 ;
      const int hmmask  = halfmod - 1 ;
      bool plan_warned = false ;
      for (int i=0; i<STAGE && covered < goal; i++) {
         ll tab[6][256] ;
         int ng = build_fast_tab(algs[i].perm, tab) ;
         fill(next.begin(), next.end(), planmask::none()) ;
         fill(nextcost.begin(), nextcost.end(), LLONG_MAX) ;
         auto update = [&](int j, int d) {
            planmask newplan = plans[j].set(i) ;
            ll nc = curcost[j] + algs[i].cost ;
            if (newplan.is_none()) {
               if (!plan_warned && plans[d].is_none()) {
                  lock_guard<mutex> lk(io_mutex) ;
                  cerr << "Warning: planmask overflow in " << key
                       << " alg " << i << "; greedy plan may be incomplete\n" ;
                  plan_warned = true ;
               }
               return ;
            }
            if (nc < nextcost[d]) { nextcost[d] = nc ; next[d] = newplan ; }
         } ;
         if (base == 4) {
            for (int j=0; j<goal; j++)
               if (!plans[j].is_none()) {
                  uint64_t v = pt->unpack[j] ; ll k = 0 ;
                  for (int g=0; g<ng; g++) { k += tab[g][v & 255] ; v >>= 8 ; }
                  update(j, pt->packhi[k >> hmshift] + pt->packlo[k & hmmask]) ;
               }
         } else if (base == 3) {
            for (int j=0; j<goal; j++)
               if (!plans[j].is_none()) {
                  ll v = pt->unpack[j] ; ll k = 0 ;
                  for (int g=0; g<ng; g++) { k += tab[g][v % 81] ; v /= 81 ; }
                  update(j, pt->packhi[k/halfmod] + pt->packlo[k%halfmod]) ;
               }
         } else { // base == 2
            for (int j=0; j<goal; j++)
               if (!plans[j].is_none()) {
                  uint64_t v = pt->unpack[j] ; ll k = 0 ;
                  for (int g=0; g<ng; g++) { k += tab[g][v & 15] ; v >>= 4 ; }
                  update(j, pt->packhi[k >> hmshift] + pt->packlo[k & hmmask]) ;
               }
         }
         for (int j=0; j<goal; j++)
            if (!next[j].is_none() && plans[j].is_none()) {
               plans[j] = next[j] ;
               curcost[j] = nextcost[j] ;
               covered++ ;
            }
      }
      finish_plan() ;
      if (goal > 100000) {
         lock_guard<mutex> lk(io_mutex) ;
         cout << "Planned greedy " << key << " in " << (walltime() - t0) << endl ;
      }
   }
   void plan_optimal() {
      ensure_algs_loaded() ;
      double t0 = walltime() ;
      // Same BFS structure as plan_greedy() so plans[j].test(i) is never true
      // during round i's inner loop (plans are only committed after the loop).
      // Differs from greedy only in the commit step:
      //   Greedy: only commits states not yet set (first coverage wins).
      //   Optimal: updates already-set plans if a cheaper plan is found.
      if (pt->unpack.empty())
         error("! unpack should not be empty") ;
      vector<planmask> next(goal, planmask::none()) ;
      vector<ll> nextcost(goal, LLONG_MAX) ;
      vector<ll> curcost(goal, LLONG_MAX) ;
      plans.resize(goal, planmask::none()) ;
      int starti = starti_for_plan() ;
      plans[starti] = planmask{} ;
      curcost[starti] = 0 ;
      const bool pow2base = !(base & (base - 1)) ;
      const int hmshift = pow2base ? __builtin_ctz(halfmod) : 0 ;
      const int hmmask  = halfmod - 1 ;
      bool plan_warned = false ;
      for (int i=0; i<STAGE; i++) {
         ll tab[6][256] ;
         int ng = build_fast_tab(algs[i].perm, tab) ;
         fill(next.begin(), next.end(), planmask::none()) ;
         fill(nextcost.begin(), nextcost.end(), LLONG_MAX) ;
         auto update = [&](int j, int d) {
            planmask newplan = plans[j].set(i) ;
            ll nc = curcost[j] + algs[i].cost ;
            if (newplan.is_none()) {
               if (!plan_warned && nc < curcost[d]) {
                  lock_guard<mutex> lk(io_mutex) ;
                  cerr << "Warning: planmask overflow in " << key
                       << " alg " << i << "; optimal plan may be suboptimal\n" ;
                  plan_warned = true ;
               }
               return ;
            }
            if (nc < nextcost[d]) { nextcost[d] = nc ; next[d] = newplan ; }
         } ;
         if (base == 4) {
            for (int j=0; j<goal; j++)
               if (!plans[j].is_none()) {
                  uint64_t v = pt->unpack[j] ; ll k = 0 ;
                  for (int g=0; g<ng; g++) { k += tab[g][v & 255] ; v >>= 8 ; }
                  update(j, pt->packhi[k >> hmshift] + pt->packlo[k & hmmask]) ;
               }
         } else if (base == 3) {
            for (int j=0; j<goal; j++)
               if (!plans[j].is_none()) {
                  ll v = pt->unpack[j] ; ll k = 0 ;
                  for (int g=0; g<ng; g++) { k += tab[g][v % 81] ; v /= 81 ; }
                  update(j, pt->packhi[k/halfmod] + pt->packlo[k%halfmod]) ;
               }
         } else { // base == 2
            for (int j=0; j<goal; j++)
               if (!plans[j].is_none()) {
                  uint64_t v = pt->unpack[j] ; ll k = 0 ;
                  for (int g=0; g<ng; g++) { k += tab[g][v & 15] ; v >>= 4 ; }
                  update(j, pt->packhi[k >> hmshift] + pt->packlo[k & hmmask]) ;
               }
         }
         for (int j=0; j<goal; j++)
            if (!next[j].is_none() && nextcost[j] < curcost[j]) {
               plans[j] = next[j] ;
               curcost[j] = nextcost[j] ;
            }
      }
      finish_plan() ;
      if (goal > 100000) {
         lock_guard<mutex> lk(io_mutex) ;
         cout << "Planned optimal " << key << " in " << (walltime() - t0) << endl ;
      }
   }
   void plan() {
      if (plans.size())
         return ;
      if (optimal)
         plan_optimal() ;
      else
         plan_greedy() ;
   }
   planmask getplan(xorbit &orb) {
      if (plans.empty())
         plan() ;
      ll k = 0 ;
      for (int i=23; i>=0; i--) {
         if (tolowdig[i] >= 0) {
            k *= base ;
            if (key[orb.a[i]] == 'x')
               error("! oops; hit an x") ;
            k += key[orb.a[i]] - '0' ;
         }
      }
      return plans[pt->packhi[k/halfmod]+pt->packlo[k%halfmod]] ;
   }
   void doplan(planmask plan, xorbit &orb) {
      if (plan.is_bitmask()) {
         for (int pi=STAGE-1; pi>=0; pi--)
            if (plan.test(pi))
               permimul(orb, algs[pi].perm) ;
      } else {
         int fields[7], n = plan.decode(fields) ;
         for (int k = n-1; k >= 0; k--)
            permimul(orb, algs[fields[k]].perm) ;
      }
   }
   FILE *openprobcache(const char *how) {
      int len = filename.size() ;
      if (len < 4 || strcmp(filename.c_str()+len-4, ".alg") != 0) {
         cerr << "Looking at filename " << filename << endl ;
         error("! algorithm filename seems bogus") ;
      }
      string probname(filename) ;
      probname.resize(len-3) ;
      probname += optimal ? "prob3-optimal" : "prob3-greedy" ;
      FILE *r = fopen(probname.c_str(), how) ;
      return r ;
   }
   // Returns true if probabilities were computed fresh (cache was missing),
   // false if loaded from an existing cache file.
   bool ensureprobs() {
      if (goal > MINWRITE) {
         FILE *f = openprobcache("r") ;
         if (f) {
            ll cost ; double v ;
            while (true) {
               int r = fscanf(f, "%lld", &cost) ;
               if (r != 1) break ;
               int r2 = fscanf(f, " %lg", &v) ;
               if (r2 == 1) {
                  algcosts.push_back(cost) ;
                  stageprob.push_back(v) ;
               } else {
                  max_activations_seen = (int)cost ;
                  break ;
               }
            }
            fclose(f) ;
            if (stageprob.size() != (size_t)STAGE) {
               cerr << "For " << filename << endl ;
               error("! probability cache file seems corrupt") ;
            }
            used_perms = 0 ;
            for (auto v : stageprob)
               if (v > 0) used_perms++ ;
         }
      }
      bool computed = (stageprob.size() == 0) ;
      if (computed) {
         plan() ;
         if (bad)
            return false ;
         if (goal > MINWRITE) {
            FILE *f = openprobcache("w") ;
            if (f == 0)
               error("! cannot write file probability cache") ;
            for (int i = 0; i < (int)stageprob.size(); i++)
               fprintf(f, "%lld %.15lg\n", algcosts[i], stageprob[i]) ;
            fprintf(f, "%d\n", max_activations_seen) ;
            fclose(f) ;
         }
      }
      plans.clear() ;
      plans.shrink_to_fit() ;
      return computed ;
   }
   reductionwplan cost(cubesizetype Mx, cubesizetype My, char typ) {
      reductionwplan r ;
      r.cost = 0 ;
      r.net = this ;
      if (typ == 'a' && stageprob.size() == 0)
         ensureprobs() ;
      if (algcosts.empty())
         ensure_algs_loaded() ;
      for (int i=0; i<(int)algcosts.size(); i++) {
         wplan rr ;
         if (typ == 'a')
            rr = ::avcost(algcosts[i], Mx, My, stageprob[i]) ;
         else if (typ == 'w')
            rr = ::wccost(algcosts[i], Mx, My) ;
         else
            error("! bad cost type") ;
         r.plan.push_back(rr) ;
         r.cost += rr.cost ;
      }
      return r ;
   }
} ;
/*
 *   Take a string indicating a partitioning or subset of the six
 *   faces (like 0x0xx0 or 012344) and reorder it in the one of the
 *   24 rotations that gives the lowest lexicographical version.
 *
 *   If mod is 1, then all numbers get smashed to zero (this lets us
 *   take a partitioning instruction like x0x11x and turn it into
 *   a face subset (like 00x0xx).
 */
string canon(string key, int mod=0) {
   string r, k2(6, 0) ;
   for (int m=0; m<24; m++) {
      int renum[6] ;
      int at = '0' ;
      for (int i=0; i<6; i++)
         renum[i] = 0 ;
      for (int i=0; i<6; i++) {
         int c = key[rots[m][i]] ;
         if (c == 'x') {
            k2[i] = c ;
         } else {
            if (renum[c-'0'] == 0) {
               renum[c-'0'] = at ;
               if (!mod)
                  at++ ;
            }
            k2[i] = renum[c-'0'] ;
         }
      }
      if (r.size() == 0 || k2 < r)
         r = k2 ;
   }
   return r ;
}
map<string, vector<netgen *>> algset ;
void addgen(string key, netgen *ng) {
   key = canon(key, 1) ;
   algset[key].push_back(ng) ;
}
/*
 *   Read a set of algorithms from an "algs" directory.  Only read files
 *   with an ".alg" suffix and whose first six characters of the filename
 *   looks like a valid key followed by a hyphen.
 */
void readalgs() {
   expandrotations() ;
   duration() ;
   const string &dirname = algsdir ;
   DIR *dir = opendir(dirname.c_str()) ;
   if (dir == 0)
      error("! no algs directory to read") ;
   int nalgsloaded = 0 ;
   while (struct dirent *entry = readdir(dir)) {
      int len = strlen(entry->d_name) ;
      if (len < 4 || strcmp(entry->d_name + len - 4, ".alg") != 0)
         continue ; // not an .alg file (e.g. .prob3-greedy, .eval, .) — silently skip
      // Everything below is an .alg file; any format problem is fatal.
      int validalgfile = 1 ;
      for (int i=0; validalgfile && i<6; i++)
         if (strchr("x0123", entry->d_name[i]) == 0)
            validalgfile = 0 ;
      if (!validalgfile || entry->d_name[6] != '-') {
         fprintf(stderr, "! invalid .alg filename (bad key or missing '-'): %s\n",
                 entry->d_name) ;
         error("! invalid .alg filename") ;
      }
      // Distinguish old format ({key}-{packed_cost}-{hash}.alg, 2 hyphens)
      // from 3-field format ({key}-{depth}-{f1}-{f2}-{f3}-{hash}.alg, 5 hyphens)
      // and 4-field format ({key}-{depth}-{v}-{h}-{o}-{c}-{hash}.alg, 6 hyphens).
      int nhyphens = 0 ;
      for (int i=6; i<len-4; i++)
         if (entry->d_name[i] == '-') nhyphens++ ;
      if (nhyphens == 5 || nhyphens == 6) {
         // depth is always the second field (chars 7 onward until next hyphen)
         long long depth = atoll(entry->d_name + 7) ;
         if (depth > planmask::max_stages()) {
            fprintf(stderr, "! .alg filename depth %lld exceeds max_stages (%d): %s\n",
                    depth, planmask::max_stages(), entry->d_name) ;
            error("! .alg filename depth exceeds max_stages") ;
         }
      } else if (nhyphens != 2) {
         fprintf(stderr, "! unrecognised .alg filename format (%d hyphens): %s\n",
                 nhyphens, entry->d_name) ;
         error("! unrecognised .alg filename format") ;
      }
      const char *key = strndup(entry->d_name, 6) ;
      string filepath = dirname + "/" + entry->d_name ;
      FILE *f = fopen(filepath.c_str(), "r") ;
      if (f == 0)
         error("! cannot read algorithm file") ;
      vector<string> algs ;
      char inpline[1024] ;
      while (fgets(inpline, 1000, f)) {
         char *p = inpline + strlen(inpline) - 1 ;
         while (p >= inpline && *p <= ' ')
            *p-- = 0 ;
         if (p > inpline)
            algs.push_back(inpline) ;
      }
      fclose(f) ;
      netgen *gen_greedy = new netgen(key, algs, filepath, 0, false) ;
      addgen(key, gen_greedy) ;
      if (!greedyonly) {
         netgen *gen_optimal = new netgen(key, algs, filepath, 0, true) ;
         addgen(key, gen_optimal) ;
      }
      nalgsloaded++ ;
   }
   closedir(dir) ;
   cout << "Loaded " << nalgsloaded << " alg" << (nalgsloaded == 1 ? "" : "s")
        << " in " << duration() << endl ;
}
map<string, double> bestcostcache ;
map<string, reductionwplan> bestcostplan ;
int splitby(const string &key, string &k2, string &k3) {
   int firstmod = -1 ;
   for (int i=0; firstmod < 0 && i<5; i++)
      if (key[i] != 'x')
         for (int j=i+1; firstmod < 0 && j<6; j++)
            if (key[j] == key[i])
               firstmod = i ;
   if (firstmod == -1)
      return 0 ; // all done
   k2 = string(6, 0) ;
   k3 = key ;
   for (int i=0; i<6; i++)
      if (key[i] == key[firstmod]) {
         k2[i] = '0' ;
         k3[i] = 'x' ;
      } else
         k2[i] = 'x' ;
   remapper remap ;
   for (int i=0; i<6; i++)
      if (k3[i] != 'x')
         k3[i] = remap.get(k3[i]) ;
   return 1 ;
}
double bestcost(string key, cubesizetype Mx, cubesizetype My, char typ) {
   // set any singletons to x before canonicalization
   int cnt[6] ;
   for (int i=0; i<6; i++)
      cnt[i] = 0 ;
   for (int i=0; i<6; i++)
      if (key[i] != 'x')
         cnt[key[i]-'0']++ ;
   for (int i=0; i<6; i++)
      if (key[i] != 'x' && cnt[key[i]-'0'] == 1)
         key[i] = 'x' ;
   auto ckey = canon(key) ;
   if (bestcostcache.find(ckey) != bestcostcache.end())
      return bestcostcache[ckey] ;
   string k2, k3 ;
   if (!splitby(key, k2, k3))
      return 0 ;
   auto remcost = bestcost(k3, Mx, My, typ) ;
   auto ck2 = canon(k2) ;
   if (algset.find(ck2) == algset.end())
      error("! can't find any good continuations") ;
   double r = INFCOST ;
   reductionwplan decomp ;
   for (auto netgp: algset[ck2]) {
      if (netgp->bad) continue ;
      auto t = netgp->cost(Mx, My, typ) ;
      auto rem2 = bestcost(netgp->key, Mx, My, typ) + remcost ;
      if (t.cost + rem2 < r) {
         r = t.cost + rem2 ;
         decomp = t ;
      }
   }
   bestcostcache[ckey] = r ;
   bestcostplan[ckey] = decomp ;
   return r ;
}
vector<reductionwplan> reductionplan ;
void explain(string key) {
   while (1) {
      string k2, k3 ;
      if (key == "012345") {
         cout << "Finished with key " << key << endl ;
         return ;
      }
      if (!splitby(key, k2, k3)) {
         cerr << "Could not spit " << key << endl ;
         error("! could not split key") ;
      }
      auto ckey = canon(k2) ;
      if (bestcostplan.find(ckey) == bestcostplan.end())
         error("! internal error; no progression?") ;
      auto plan = &bestcostplan[ckey] ;
      string newkey(6, 0) ;
      int found = 0 ;
      /*
       *   Here we need to search for a rotation that lets us apply a given
       *   netgen in a prescribed way.  We have the original partitioning,
       *   the destination partitioning, and the netgen; try to find the
       *   appropriate steps.
       */
      int m ;
      for (m=0; m<24; m++) {
         remapper remap ;
         for (int i=0; i<6; i++)
            newkey[rots[m][i]] = remap.get(plan->net->key[i]) ;
         int good = 1 ;
         int remapt[6] ;
         for (int i=0; i<6; i++)
            remapt[i] = -1 ;
         for (int i=0; i<6; i++)
            if (k2[i] != 'x') {
               if (newkey[i] == 'x')
                  good = 0 ;
               else {
                  if (remapt[k2[i]-'0'] < 0)
                     remapt[k2[i]-'0'] = k2[i] ;
                  else if (remapt[k2[i]-'0'] != k2[i])
                     good = 0 ;
               }
            } else {
               if (newkey[i] != 'x')
                  good = 0 ;
            }
         if (good) {
            found = 1 ;
            break ;
         }
      }
      if (!found) {
         cerr << "Trying to use " << plan->net->key << " for " << key << endl ;
         error("! could not find right rotation") ;
      }
      const char *mode = plan->net->optimal ? "optimal" : "greedy" ;
      if (plan->net->used_perms < 0)
         plan->net->plan() ;
      cout << "Use " << newkey << " " << plan->net->filename << " [" << mode << "] cost " << plan->cost << " on " << key << " map " << m << " used " << plan->net->used_perms << "/" << plan->net->stringalgs.size() << endl ;
      netgen *f = plan->net ;
      if (m != 0)
         f = new netgen(f->key, f->stringalgs, f->filename, m, f->optimal) ;
      reductionplan.push_back({plan->cost, f, plan->plan}) ;
      string k4(6, 0) ;
      remapper remap ;
      for (int i=0; i<6; i++)
         if (newkey[i] == 'x')
            k4[i] = remap.get(key[i]) ;
         else
            k4[i] = remap.get(6+newkey[i]) ;
      swap(key, k4) ;
   }
}
void build_wing_cycle_index() {
   for (int a=0; a<24; a++)
      for (int b=0; b<24; b++)
         for (int c=0; c<24; c++)
            wing_cycle_index[a][b][c] = -1 ;
   int idx = 0 ;
   for (int a=0; a<24; a++)
      for (int b=0; b<24; b++)
         for (int c=0; c<24; c++) {
            if (a==b || b==c || a==c)
               continue ;
            if (a < b && a < c) {            // canonical
               wing_cycle_index[a][b][c] = idx ;
               wing_cycle_index[b][c][a] = idx ; // same cycle, other rotation
               wing_cycle_index[c][a][b] = idx ; // same cycle, other rotation
               wing_cycle_triple[idx][0] = a ;
               wing_cycle_triple[idx][1] = b ;
               wing_cycle_triple[idx][2] = c ;
               idx++ ;
            }
         }
   if (idx != 4048)
      error("error, not 4048") ;
}

vector<int> solve_wings(vector<int> orb) {
   vector<int> algs ;
   for (int i=0; i<24; i++) {
      while (orb[i] != i) { // this skips solved slots
         int t1 = find(orb.begin(), orb.end(), i) - orb.begin(); // slot where piece i currently is
         int t2 = i ; // next target is the spot where the first unsolved piece should go
         int t3 = orb[i] ; // the target of the piece at t2
         if (t3 == t1) { // to resolve a two-cycle, choose a victim piece that is part of an even cycle.
            bool seen[24] = {} ;
            int cyclebreak = -1 ;
            for (int s=i+1; s<24 && cyclebreak < 0; s++) {
               if (orb[s] == s || s == t1 || seen[s]) 
                  continue ;
               int len = 0;
               int x = s ;
               do { seen[x] = true ; len++ ; x = orb[x] ; } while (x != s) ; // count the length of the cycle
               if ((len % 2) == 0) {
                  cyclebreak = s ;
               }
            }
            t3 = cyclebreak ;
         }
         algs.push_back(wing_cycle_index[t1][t2][t3]) ;
         int tmp = orb[t3] ;
         orb[t3] = orb[t2] ;
         orb[t2] = orb[t1] ;
         orb[t1] = tmp ;
      }
   }
   return algs ;
}
vector<mv> wing_alg_moves[4048] ;
int wing_cycle_of(const vector<mv> &moves) {
   xcube cube(4, 5, 'U') ;
   for (auto mv : moves)
      domove(mv, cube) ;
   xorbit &o = cube.wingorbits[0] ;
   int dest[24] ;
   for (int t=0; t<24; t++) 
      dest[o.a[t]] = t ;
   int t1 = -1;
   int moved = 0 ;
   for (int s=0; s<24; s++)
      if (dest[s] != s) { 
         moved++ ; 
         if (t1 < 0) 
            t1 = s ; 
      }
   if (moved != 3) 
      return -1 ;
   return wing_cycle_index[t1][dest[t1]][dest[dest[t1]]] ;
}

vector<mv> rotate_moves(const vector<mv> &in, int m, bool inv) {
   vector<mv> r = in ;
   if (inv) {
      reverse(r.begin(), r.end()) ;
      for (auto &mv : r) 
         mv.twist = (4 - mv.twist) & 3 ;
   }
   for (auto &mv : r) 
      mv.face = rots[m][mv.face] ;
   return r ;
}
const char *wing_base_algs[] = {
   "F D 2B' D' F' D 2B D'",
   "F' L' 2B L F L' 2B' L",
   "F D 2B2 D' F' D 2B2 D'",
   "F' L' 2B2 L F L' 2B2 L",
   "F D 2B D' F' D 2B' D'",
   "F' L' 2B' L F L' 2B L",
   "F D 2F' D' F' D 2F D'",
   "F' L' 2F L F L' 2F' L",
   "F D 2F2 D' F' D 2F2 D'",
   "F' L' 2F2 L F L' 2F2 L",
   "F D2 F' 2D F D2 F' 2D'",
   "F' L2 F 2L' F' L2 F 2L",
   "F D2 F' 2D2 F D2 F' 2D2",
   "F' L2 F 2L2 F' L2 F 2L2",
   "F D2 F' 2D' F D2 F' 2D",
   "F' L2 F 2L F' L2 F 2L'",
   "F D2 F' 2U F D2 F' 2U'",
   "F' L2 F 2R' F' L2 F 2R",
   "F D' 2B2 D F' D' 2B2 D",
   "F' L 2B2 L' F L 2B2 L'",
   "F D2 F' 2U2 F D2 F' 2U2",
   "F' L2 F 2R2 F' L2 F 2R2",
   "F D' 2F2 D F' D' 2F2 D",
   "F' L 2F2 L' F L 2F2 L'",
   "F D2 F' 2U' F D2 F' 2U",
   "F' L2 F 2R F' L2 F 2R'",
   "F D' F' 2D2 F D F' 2D2",
   "F' L F 2L2 F' L' F 2L2",
   "F D' F' 2D' F D F' 2D",
   "F' L F 2L F' L' F 2L'",
   "F D' F' 2U F D F' 2U'",
   "F' L F 2R' F' L' F 2R",
   "F D' F' 2U2 F D F' 2U2",
   "F' L F 2R2 F' L' F 2R2",
   "F 2D F' D2 F 2D' F' D2",
   "F' 2L' F L2 F' 2L F L2",
   "F 2D F' U2 F 2D' F' U2",
   "F' 2L' F R2 F' 2L F R2",
   "F 2D2 F' D2 F 2D2 F' D2",
   "F' 2L2 F L2 F' 2L2 F L2",
   "F 2D2 F' U2 F 2D2 F' U2",
   "F' 2L2 F R2 F' 2L2 F R2",
   "F 2D' F' U2 F 2D F' U2",
   "F' 2L F R2 F' 2L' F R2",
   "F D 2B D' F D 2B' D' F2",
   "F' L' 2B' L F' L' 2B L F2",
   "F D 2B' D' F D 2B D' F2",
   "F D 2F D' F D 2F' D' F2",
   "F D 2B D' F2 D 2B' D' F",
   "F' L' 2B' L F2 L' 2B L F'",
   "F D F 2D F' D' F 2D' F2",
   "F' L' F' 2L' F L F' 2L F2",
   "F D2 F 2D F' D2 F 2D' F2",
   "F' L2 F' 2L' F L2 F' 2L F2",
   "F D' F 2D F' D F 2D' F2",
   "F' L F' 2L' F L' F' 2L F2",
   "F D' F 2D2 F' D F 2D2 F2",
   "F' L F' 2L2 F L' F' 2L2 F2",
   "F 2D F D2 F' 2D' F D2 F2",
   "F 2D F U F' 2D' F U' F2",
   "F' 2L' F' R' F 2L F' R F2",
   "F 2D2 F D F' 2D2 F D' F2",
   "F 2D2 F U F' 2D2 F U' F2",
   "F B D F' 2D' F D' F' 2D B'",
   "F B D F' 2U F D' F' 2U' B'",
   "F B D' F' 2D' F D F' 2D B'",
   "F' B' L F 2L F' L' F 2L' B",
   "F B 2D2 F' D F 2D2 F' D' B'",
   "F' B' 2L2 F L' F' 2L2 F L B",
   "F B' D F' 2D' F D' F' 2D B",
   "F' B L' F 2L F' L F 2L' B'",
   "F D R' 2D2 R D' R' 2D2 R F'",
   "F D R' 2U2 R D' R' 2U2 R F'",
   "F D2 L' 2D' L D2 L' 2D L F'",
   "F' L2 D 2L D' L2 D 2L' D' F",
   "F D B2 L' 2F L B2 L' 2F' L D' F'",
   "F D' L F' D 2L D' F L' D F' 2L'",
   "F' L D' F L' 2D' L F' D L' F 2D",
   "F D' L F' D 2L2 D' F L' D F' 2L2",
   "F D' L 2F2 L' D F' L D' 2F2 D L'",
   "F' L D' 2F2 D L' F D' L 2F2 L' D",
   "F D' L 2F' L' D F' L D' 2F D L'",
   "F' L D' 2F D L' F D' L 2F' L' D",
   "F D L' F 2D' F' L D' F L' 2D L F2",
   "F' L' D F' 2L F D' L F' D 2L' D' F2",
   "F D L' F 2U F' L D' F L' 2U' L F2",
   "F' L' D F' 2R' F D' L F' D 2R D' F2",
} ;
void loadwingalgs() {
   for (int idx=0; idx<4048; idx++)
      wing_alg_moves[idx].clear() ;
   for (int i=0; i<87; i++) {
      vector<mv> m0 = parse_moves(wing_base_algs[i]) ;
      for (int m=0; m<24; m++)
         for (int inv=0; inv<2; inv++) {
            vector<mv> v = rotate_moves(m0, m, inv) ;
            int idx = wing_cycle_of(v) ;
            if (idx >= 0 && wing_alg_moves[idx].empty())
               wing_alg_moves[idx] = v ;
         }
   }
   for (int idx=0; idx<4048; idx++)
      if (wing_alg_moves[idx].empty())
         error("! wing base set does not cover all 4048 cycles") ;
}
/*
 *   Execute one wing 3-cycle (cycle index idx) across a set of orbit depths at
 *   once.  The alg's depth-2 inner slice is expanded into one parallel slice per
 *   target orbit.
 */
ll doblock_wings(xcube &cube, int idx, const vector<int> &depths) {
#if WRITEFRAMES
   if (writemoveinc && runningmovecount >= writemovetarget)
      writeframe(cube) ;
#endif
   blocks++ ;
   auto &alg = wing_alg_moves[idx] ;
   int nslice = 0 ;
   for (auto &mv : alg)
      if (mv.dep == 2)
         nslice++ ;
   if (writemoves) {
      for (auto mv : alg) {
         if (mv.dep == 2) {
            for (int k : depths) { 
               mv.dep = k + 2 ; 
               curblock.push_back(mv) ; 
            }
         } else {
            curblock.push_back(mv) ;
         }
      }
      sendblock() ;
   }
   for (auto mv : alg) {
      if (mv.dep == 2) {
         for (int k : depths) { 
            mv.dep = k + 2 ; 
            domove(mv, cube) ; 
         }
      } else {
         domove(mv, cube) ;
      }
   }
   ll mvs = (ll)alg.size() + (ll)nslice * ((ll)depths.size() - 1) ;
   runningmovecount += mvs ;
   return mvs ;
}
void solvecorners(xcube &cube) {
   int odd = (int)((cube.Mx + cube.My) & 1) ;
   threecube cc = cube.tc ;
   static int inited = 0 ;
   vector<mv> cornersolution ;
   if (odd) {
      ll midslice = (ll)(cube.Mx + 1) ;
      // orient center centers first.
      for (int mvs=0; mvs<2; mvs++) {
         int centerwrongcount = cc.countbadcenters() ;
         if (centerwrongcount) {
            int found = 0 ;
            for (int f=0; found==0 && f<3; f++)
               for (int tw=1; found==0 && tw<4; tw++) {
                  cc.domove({2, (short)f, (short)tw}) ;
                  if (cc.countbadcenters() < centerwrongcount) {
                     found = 1 ;
                     cornersolution.push_back({(int)midslice, (short)f, (short)tw}) ;
                  } else {
                     cc.domove({2, (short)f, (short)(4-tw)}) ;
                  }
               }
         }
      }
      if (cc.countbadcenters())
         error("! impossible centercenter state") ;
      if (cornersolution.size()) {
         cout << "Oriented center centers in " << cornersolution.size() << " moves." << endl ;
         for (auto &mv: cornersolution)
            domove(mv, cube) ;
      }
      static struct prunetable<7, 10000019, 12> phase1prune("phase1") ;
      static struct prunetable<8, 50000017, 18, 1> phase2prune("phase2") ;
      if (!inited) {
         inited = 1 ;
         threecube c1('1', 1) ;
         phase1prune.fillit(c1) ;
         threecube c2('C', 1) ;
         phase2prune.fillit(c2) ;
      }
      threecube c1 = cc.tophase1() ;
      int r1 = phase1prune.solve(c1) ;
      if (r1 < 0)
         error("! failed to solve phase 1") ;
      cornersolution.insert(cornersolution.end(),
                            phase1prune.mvs.begin(), phase1prune.mvs.end()) ;
      for (auto &mv: phase1prune.mvs)
         domove(mv, cube) ;
      cc = cube.tc ;
      int r2 = phase2prune.solve(cc) ;
      if (r2 < 0)
         error("! failed to solve phase 2") ;
      for (auto &mv: phase2prune.mvs)
         domove(mv, cube) ;
      cornersolution.insert(cornersolution.end(),
                            phase2prune.mvs.begin(), phase2prune.mvs.end()) ;
      cout << "Solved corners plus middle edges in " << r1 << " + " << r2
                               << " = " << (r1 + r2) << " moves." << endl ;
   } else {
      static struct prunetable<6, 2000003, 11> cornerprune("corner") ;
      if (!inited) {
         inited = 1 ;
         threecube c1('C', 0) ;
         cornerprune.fillit(c1) ;
      }
      int r = cornerprune.solve(cc) ;
      if (r < 0)
         error("! failed to solve corners") ;
      for (auto &mv: cornerprune.mvs)
         domove(mv, cube) ;
      cornersolution = cornerprune.mvs ;
      cout << "Solved corners in " << cornerprune.mvs.size() << " moves." << endl ;
   }
   runningmovecount += cornersolution.size() ;
   if (writemoves) {
      for (auto mv: cornersolution) {
         cout << " " ;
         showmove(mv) ;
      }
      cout << endl ;
   }
}
// Solve every wing orbit with 3-cycles, executed in wing-cycle-index order
ll do_wing_cycles(xcube &cube) {
   ll wr = 0 ;
   vector<vector<int>> seq((size_t)cube.My) ;
   for (ll k=0; k<cube.My; k++)
      seq[(size_t)k] = solve_wings(cube.get_wings((ll)k)[0]) ;
   // Index-ordered execution.
   vector<vector<int>> bucket(4048) ;
   for (ll k=0; k<cube.My; k++)
      for (int idx : seq[k])
         bucket[idx].push_back(k) ;
   for (int v=0; v<4048; v++)
      if (!bucket[v].empty()) {
         ll wm = doblock_wings(cube, v, bucket[v]) ;
         wr += wm ;
         wingmoves += wm ;
      }
   return wr ;
}
// Strip-local rotation helpers for fill_strip_from_moveset.
// srows[0..pw-1][0..Mx-1]; fr[6] is current accumulated face rotation.
static void s_rot4hor(vector<vector<xorbit>> &srows, int f1, int o1, int f2, int o2,
                      int f3, int o3, int f4, int o4, int local_dep, const int fr[6]) {
   for (auto &orb : srows[local_dep]) {
      uchar &c1 = orb.a[4*f1+((o1+fr[f1])&3)] ;
      uchar &c2 = orb.a[4*f2+((o2+fr[f2])&3)] ;
      uchar &c3 = orb.a[4*f3+((o3+fr[f3])&3)] ;
      uchar &c4 = orb.a[4*f4+((o4+fr[f4])&3)] ;
      auto t = c1 ; c1 = c2 ; c2 = c3 ; c3 = c4 ; c4 = t ;
   }
}
static void s_rot2hor(vector<vector<xorbit>> &srows, int f1, int o1, int f2, int o2,
                      int f3, int o3, int f4, int o4, int local_dep, const int fr[6]) {
   for (auto &orb : srows[local_dep]) {
      uchar &c1 = orb.a[4*f1+((o1+fr[f1])&3)] ;
      uchar &c2 = orb.a[4*f2+((o2+fr[f2])&3)] ;
      uchar &c3 = orb.a[4*f3+((o3+fr[f3])&3)] ;
      uchar &c4 = orb.a[4*f4+((o4+fr[f4])&3)] ;
      swap(c1, c3) ;
      swap(c2, c4) ;
   }
}
static void s_rot4ver(vector<vector<xorbit>> &srows, int f1, int o1, int f2, int o2,
                      int f3, int o3, int f4, int o4, int dep, const int fr[6]) {
   for (auto &row : srows) {
      uchar &c1 = row[dep].a[4*f1+((o1+fr[f1])&3)] ;
      uchar &c2 = row[dep].a[4*f2+((o2+fr[f2])&3)] ;
      uchar &c3 = row[dep].a[4*f3+((o3+fr[f3])&3)] ;
      uchar &c4 = row[dep].a[4*f4+((o4+fr[f4])&3)] ;
      auto t = c1 ; c1 = c2 ; c2 = c3 ; c3 = c4 ; c4 = t ;
   }
}
static void s_rot2ver(vector<vector<xorbit>> &srows, int f1, int o1, int f2, int o2,
                      int f3, int o3, int f4, int o4, int dep, const int fr[6]) {
   for (auto &row : srows) {
      uchar &c1 = row[dep].a[4*f1+((o1+fr[f1])&3)] ;
      uchar &c2 = row[dep].a[4*f2+((o2+fr[f2])&3)] ;
      uchar &c3 = row[dep].a[4*f3+((o3+fr[f3])&3)] ;
      uchar &c4 = row[dep].a[4*f4+((o4+fr[f4])&3)] ;
      swap(c1, c3) ;
      swap(c2, c4) ;
   }
}

struct StripSolver {
   vector<mv> prevmvs ;
   ll cancellations = 0 ;
   ll block_count = 0 ;
   vector<vector<planmask>> movemasks ;
   vector<ll> bits ;
   vector<int> nxt, prv ;
   vector<int> xs, ys ;
   vector<reductionwplan> localplan ;
   mt19937_64 local_rng ;
   struct PhaseStat { ll hits = 0, moves = 0 ; double expectedcost = 0 ; } ;
   vector<PhaseStat> phase_stats ;

   StripSolver(ll phasewidth_, ll Mx_) : local_rng(random_device{}()) {
      movemasks.resize((size_t)phasewidth_) ;
      for (auto &mm : movemasks)
         mm.resize((size_t)Mx_) ;
      localplan = reductionplan ;
   }

   void fill_row(xcube &cube, ll row) {
      cube.rows[(size_t)row].resize((size_t)cube.Mx) ;
      for (auto &orb : cube.rows[(size_t)row]) {
         oclear(orb) ;
         shuffle(orb.a, orb.a+24, local_rng) ;
      }
   }

   // Apply one pre-normalized move to the strip-local srows and fr.
   // dep < My_: normal slice (hor if in strip, ver across all strip rows).
   // dep >= My_: middle slice (ver only, no ccs).
   // dep == -1 (face move): update fr only.
   static void apply_strip_move(const mv &smv, int fr[6],
                                vector<vector<xorbit>> &srows,
                                ll strip_start, ll pw, ll My_) {
      if (smv.dep == 1) {
         fr[smv.face] = (fr[smv.face] + 4 - smv.twist) & 3 ;
         return ;
      }
      int dep = smv.dep - 2 ;
      int local_dep = dep - (int)strip_start ;
#define SH(fn,a,b,c,d,e,f,g,h) if(local_dep>=0&&local_dep<(int)pw) s_##fn##hor(srows,a,b,c,d,e,f,g,h,local_dep,fr)
#define SV(fn,a,b,c,d,e,f,g,h) s_##fn##ver(srows,a,b,c,d,e,f,g,h,dep,fr)
#define SV2(fn,a,b,c,d,e,f,g,h) s_##fn##ver(srows,a,b,c,d,e,f,g,h,dep,fr);s_##fn##ver(srows,a,b^2,c,d^2,e,f^2,g,h^2,dep,fr)
      if (dep >= (int)My_) {
         switch (smv.face * 3 + smv.twist) {
case 1:  SV2(rot4,1,1,2,1,3,1,4,1); break;
case 2:  SV2(rot2,1,1,2,1,3,1,4,1); break;
case 3:  SV2(rot4,4,1,3,1,2,1,1,1); break;
case 4:  SV2(rot4,2,0,0,0,4,2,5,0); break;
case 5:  SV2(rot2,2,0,5,0,4,2,0,0); break;
case 6:  SV2(rot4,2,0,5,0,4,2,0,0); break;
case 7:  SV2(rot4,0,3,1,2,5,1,3,0); break;
case 8:  SV2(rot2,0,3,3,0,5,1,1,2); break;
case 9:  SV2(rot4,0,3,3,0,5,1,1,2); break;
case 10: SV2(rot4,2,2,5,2,4,0,0,2); break;
case 11: SV2(rot2,2,2,5,2,4,0,0,2); break;
case 12: SV2(rot4,2,2,0,2,4,0,5,2); break;
case 13: SV2(rot4,0,1,3,2,5,3,1,0); break;
case 14: SV2(rot2,0,1,3,2,5,3,1,0); break;
case 15: SV2(rot4,0,1,1,0,5,3,3,2); break;
case 16: SV2(rot4,4,3,3,3,2,3,1,3); break;
case 17: SV2(rot2,1,3,2,3,3,3,4,3); break;
case 18: SV2(rot4,4,3,1,3,2,3,3,3); break;
         }
      } else {
         switch (smv.face * 3 + smv.twist) {
case 1:  SH(rot4,1,0,2,0,3,0,4,0); SV(rot4,1,1,2,1,3,1,4,1); break;
case 2:  SH(rot2,1,0,2,0,3,0,4,0); SV(rot2,1,1,2,1,3,1,4,1); break;
case 3:  SH(rot4,4,0,3,0,2,0,1,0); SV(rot4,4,1,3,1,2,1,1,1); break;
case 4:  SH(rot4,2,3,0,3,4,1,5,3); SV(rot4,2,0,0,0,4,2,5,0); break;
case 5:  SH(rot2,2,3,5,3,4,1,0,3); SV(rot2,2,0,5,0,4,2,0,0); break;
case 6:  SH(rot4,2,3,5,3,4,1,0,3); SV(rot4,2,0,5,0,4,2,0,0); break;
case 7:  SH(rot4,0,2,1,1,5,0,3,3); SV(rot4,0,3,1,2,5,1,3,0); break;
case 8:  SH(rot2,0,2,3,3,5,0,1,1); SV(rot2,0,3,3,0,5,1,1,2); break;
case 9:  SH(rot4,0,2,3,3,5,0,1,1); SV(rot4,0,3,3,0,5,1,1,2); break;
case 10: SH(rot4,2,1,5,1,4,3,0,1); SV(rot4,2,2,5,2,4,0,0,2); break;
case 11: SH(rot2,2,1,5,1,4,3,0,1); SV(rot2,2,2,5,2,4,0,0,2); break;
case 12: SH(rot4,2,1,0,1,4,3,5,1); SV(rot4,2,2,0,2,4,0,5,2); break;
case 13: SH(rot4,0,0,3,1,5,2,1,3); SV(rot4,0,1,3,2,5,3,1,0); break;
case 14: SH(rot2,0,0,3,1,5,2,1,3); SV(rot2,0,1,3,2,5,3,1,0); break;
case 15: SH(rot4,0,0,1,3,5,2,3,1); SV(rot4,0,1,1,0,5,3,3,2); break;
case 16: SH(rot4,4,2,3,2,2,2,1,2); SV(rot4,4,3,3,3,2,3,1,3); break;
case 17: SH(rot2,1,2,2,2,3,2,4,2); SV(rot2,1,3,2,3,3,3,4,3); break;
case 18: SH(rot4,4,2,1,2,2,2,3,2); SV(rot4,4,3,1,3,2,3,3,3); break;
         }
      }
#undef SH
#undef SV
#undef SV2
   }

   // dosolve convention: std::rotate (LEFT rotate by fr[f])
   static void apply_fr_dosolve(vector<vector<xorbit>> &srows, int fr[6]) {
      for (auto &row : srows)
         for (auto &orb : row)
            for (int f = 0; f < 6; f++)
               if (fr[f] & 3)
                  rotate(orb.a + 4*f, orb.a + 4*f + (fr[f] & 3), orb.a + 4*f + 4) ;
      for (int f = 0; f < 6; f++) fr[f] = 0 ;
   }

   // normalizefacerots convention: rot4/rot2/rot4i (RIGHT rotate by fr[f])
   static void apply_fr_normfr(vector<vector<xorbit>> &srows, int fr[6]) {
      for (auto &row : srows)
         for (auto &orb : row)
            for (int f = 0; f < 6; f++) {
               switch (fr[f] & 3) {
               case 0: break ;
               case 1: orb.rot4(4*f, 4*f+1, 4*f+2, 4*f+3) ; break ;
               case 2: orb.rot2(4*f, 4*f+1, 4*f+2, 4*f+3) ; break ;
               case 3: orb.rot4i(4*f, 4*f+1, 4*f+2, 4*f+3) ; break ;
               }
            }
      for (int f = 0; f < 6; f++) fr[f] = 0 ;
   }

   // Replay pre-normalized movesets into thread-local srows, then move to cube.rows.
   // Uses separate fr arrays per phase to match the serial path's two-step convention:
   //   - orig_mvs + dosolve: std::rotate(fr_orig)
   //   - prep_before + normalizefacerots: rot4/rot2/rot4i(fr_corner)
   //   - prep_after: wing cycles (zero net face rotation, so no final apply needed)
   void fill_strip_from_moveset(xcube &cube, ll strip_start, ll pw,
                                const vector<mv> &orig_mvs,
                                const vector<mv> &prep_before,
                                const vector<mv> &prep_after) {
      vector<vector<xorbit>> srows((size_t)pw,
                                   vector<xorbit>((size_t)cube.Mx)) ;
      for (auto &row : srows)
         for (auto &orb : row)
            oclear(orb) ;
      // Phase 1: original moveset (dosolve convention)
      int fr_orig[6] = {0,0,0,0,0,0} ;
      for (const auto &smv : orig_mvs)
         apply_strip_move(smv, fr_orig, srows, strip_start, pw, cube.My) ;
      apply_fr_dosolve(srows, fr_orig) ;  // matches dosolve's std::rotate
      // Phase 2: corner-fix moves (normalizefacerots convention)
      int fr_corner[6] = {0,0,0,0,0,0} ;
      for (const auto &smv : prep_before)
         apply_strip_move(smv, fr_corner, srows, strip_start, pw, cube.My) ;
      apply_fr_normfr(srows, fr_corner) ;  // matches normalizefacerots's rot4/rot2/rot4i
      // Phase 3: parity fix + wing cycles (net zero face rotation, no final apply)
      int fr_after[6] = {0,0,0,0,0,0} ;
      for (const auto &smv : prep_after)
         apply_strip_move(smv, fr_after, srows, strip_start, pw, cube.My) ;
      // fr_after should be 0 here (wing algos are balanced); apply anyway to be safe
      apply_fr_dosolve(srows, fr_after) ;
      for (ll wi = 0; wi < pw; wi++)
         cube.rows[(size_t)(strip_start + wi)] = std::move(srows[(size_t)wi]) ;
   }

   void release_row(xcube &cube, ll row) {
      cube.checkrow((int)row) ;
      vector<xorbit> empty ;
      swap(empty, cube.rows[(size_t)row]) ;
   }

   void trycancellation_local(mv &m) {
      if (prevmvs.empty()) {
         prevmvs.push_back(m) ;
         return ;
      }
      if (m.dep > 1) {
         prevmvs.clear() ;
         return ;
      }
      mv merged ;
      int r = trymerge(prevmvs.back(), m, merged) ;
      if (r == 2)
         prevmvs.push_back(m) ;
      else if (r == 1) {
         cancellations++ ;
         prevmvs.back() = merged ;
      } else {
         cancellations += 2 ;
         prevmvs.pop_back() ;
      }
   }

   ll doblock_fast(blockalgo &alg) {
      ll oldcancellations = cancellations ;
      block_count++ ;
      for (int i = 0; i < alg.c; i++)
         trycancellation_local(alg.mvs[i]) ;
      prevmvs.clear() ;
      for (int i = alg.c - 1; i >= 0; i--)
         trycancellation_local(alg.mvs[alg.mvs.size() - 1 - i]) ;
      return alg.o + alg.v * (ll)xs.size() + alg.h * (ll)ys.size() + alg.c -
             (cancellations - oldcancellations) ;
   }

   ll solve_strip(xcube &cube, ll strip_start, ll pw,
                  const vector<mv> *norm_mvs = nullptr,
                  const vector<mv> *prep_before = nullptr,
                  const vector<mv> *prep_after = nullptr) {
      ll r = 0 ;
      if (!alwaysinmem) {
         if (norm_mvs)
            fill_strip_from_moveset(cube, strip_start, pw,
                                    *norm_mvs, *prep_before, *prep_after) ;
         else
            for (ll wi = 0; wi < pw; wi++)
               fill_row(cube, strip_start + wi) ;
      }
      int phasebase = 0 ;
      for (auto &f : localplan) {
         for (ll ii = 0; ii < pw; ii++)
            for (ll j = 0; j < (ll)cube.Mx; j++) {
               movemasks[(size_t)ii][(size_t)j] = f.net->getplan(cube.rows[strip_start + ii][j]) ;
               f.net->doplan(movemasks[(size_t)ii][(size_t)j], cube.rows[strip_start + ii][j]) ;
            }
         for (int pi = f.net->STAGE - 1; pi >= 0; pi--, phasebase++) {
            auto &alg = f.net->algs[pi] ;
            auto &wplan = f.plan[pi] ;
            wplan.plan((cubesizetype)pw) ;
            if (stats) {
               if (phasebase >= (int)phase_stats.size())
                  phase_stats.resize(phasebase + 1) ;
               phase_stats[phasebase].expectedcost += wplan.cost ;
            }
            ll ii = 0 ;
            for (ll iii = 0; iii < (ll)(wplan.count1 + wplan.count2); iii++) {
               int w = iii < wplan.count1 ? wplan.w1 : wplan.w2 ;
               bits.assign((size_t)cube.Mx, 0) ;
               for (int wi = 0; wi < w; wi++)
                  for (ll j = 0; j < (ll)cube.Mx; j++)
                     if (movemasks[(size_t)(ii + wi)][(size_t)j].test(pi))
                        bits[(size_t)j] |= 1LL << wi ;
               prv.assign(1LL << w, -1) ;
               nxt.resize((size_t)cube.Mx) ;
               for (ll i = 0; i < (ll)cube.Mx; i++) {
                  nxt[(size_t)i] = prv[(size_t)bits[(size_t)i]] ;
                  prv[(size_t)bits[(size_t)i]] = (int)i ;
               }
               for (ll bp = 1; bp < (1LL << w); bp++) {
                  if (prv[bp] < 0)
                     continue ;
                  ys.clear() ;
                  for (int wi = 0; wi < w; wi++)
                     if ((bp >> wi) & 1)
                        ys.push_back((int)(strip_start + ii + wi)) ;
                  xs.clear() ;
                  for (int i = prv[bp]; i >= 0; i = nxt[i])
                     xs.push_back(i) ;
                  ll mvs = doblock_fast(alg) ;
                  if (stats) {
                     phase_stats[phasebase].hits += xs.size() * ys.size() ;
                     phase_stats[phasebase].moves += mvs ;
                  }
                  r += mvs ;
               }
               ii += w ;
            }
         }
      }
      for (ll wi = 0; wi < pw; wi++)
         cube.checkrelease(strip_start + wi) ;
      return r ;
   }
} ;
ll solve(xcube &cube, cubesizetype phasewidth) {
   if (phasewidth > cube.My)
      phasewidth = cube.My ;
   ll r = 0 ;
   // For streaming parallel solve, capture all prep moves so each strip can replay them.
   const bool streaming = (!alwaysinmem && readmoves) ;
   vector<mv> prep_before_norm, prep_after_norm ;
   // g_prep_recorder is always cleared to nullptr (line below) before these locals go out of scope.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdangling-pointer"
#endif
   if (streaming) g_prep_recorder = &prep_before_norm ;
   solvecorners(cube) ;
   cube.normalizefacerots() ;
   if (streaming) g_prep_recorder = &prep_after_norm ;
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
   wingmoves = 0 ;
   bool dowings = !extrapolate ; // previously had alwaysinmem condition
   // Wing parity fix
   if (dowings){
      for (ll k=0; k<cube.My; k++) {
         xorbit p = cube.wingorbits[k] ;
         int par = 0 ;
         for (int i=0; i<24; i++)
            while (p.a[i] != i) {
               swap(p.a[i], p.a[p.a[i]]) ; // this flips through the array to check its parity
               par = !par ;
            }
         if (par) {
            mv pm = { (int)(k+2), 0, 1 } ;       // quarter-turn U slice at depth k+2
            if (writemoves)
               curblock.push_back(pm) ;
            domove(pm, cube) ;
            r++ ;
            wingmoves++ ;
            runningmovecount++ ;
         }
      }
   }
   if (dowings){
      r += do_wing_cycles(cube) ;
   }
   if (streaming) g_prep_recorder = nullptr ;
   if (fast && nprob_threads > 1) {
      for (auto &f : reductionplan)
         f.net->plan() ;
      // Pre-normalize moveset for streaming replay (readmoves + parallel only)
      vector<mv> norm_moveset ;
      if (!alwaysinmem && readmoves) {
         norm_moveset = moveset ;
         for (auto &smv : norm_moveset) {
            if (smv.dep > cube.Mx + cube.My + 3 - smv.dep) {
               smv.dep = cube.Mx + cube.My + 3 - smv.dep ;
               smv.face = oppface[smv.face] ;
               smv.twist = (4 - smv.twist) & 3 ;
            }
         }
      }
      ll num_strips = ((ll)cube.My + (ll)phasewidth - 1) / (ll)phasewidth ;
      atomic<ll> next_si{0} ;
      atomic<ll> total_mvs{0} ;
      atomic<ll> total_cancellations{0} ;
      atomic<ll> total_blocks{0} ;
      if (stats) {
         int phasebase = 0 ;
         for (auto &f : reductionplan) {
            for (int pi = f.net->STAGE-1; pi >= 0; pi--, phasebase++) {
               auto &wp = f.plan[pi] ;
               wp.plan((cubesizetype)phasewidth) ;
               if (phasebase >= (int)byphase.size()) {
                  byphase.resize(phasebase+1) ;
                  byphase[phasebase].w1 = wp.w1 ;
                  byphase[phasebase].w2 = wp.w2 ;
                  byphase[phasebase].count1 = wp.count1 ;
                  byphase[phasebase].count2 = wp.count2 ;
               }
               byphase[phasebase].ba = &(f.net->algs[pi]) ;
            }
         }
      }
      // non-null only for readmoves streaming; random mode uses fill_row (null)
      const vector<mv> *norm_mvs_ptr  = streaming ? &norm_moveset    : nullptr ;
      const vector<mv> *prep_bef_ptr  = streaming ? &prep_before_norm : nullptr ;
      const vector<mv> *prep_aft_ptr  = streaming ? &prep_after_norm  : nullptr ;
      mutex stats_mu ;
      auto worker = [&]() {
         StripSolver ss((ll)phasewidth, cube.Mx) ;
         for (;;) {
            ll si = next_si.fetch_add(1) ;
            if (si >= num_strips) break ;
            ll strip_start = si * (ll)phasewidth ;
            ll pw = (ll)phasewidth ;
            if (strip_start + pw > (ll)cube.My)
               pw = (ll)cube.My - strip_start ;
            total_mvs.fetch_add(ss.solve_strip(cube, strip_start, pw,
                                               norm_mvs_ptr, prep_bef_ptr, prep_aft_ptr)) ;
         }
         total_cancellations.fetch_add(ss.cancellations) ;
         total_blocks.fetch_add(ss.block_count) ;
         if (stats) {
            lock_guard<mutex> lk(stats_mu) ;
            for (int pb = 0; pb < (int)ss.phase_stats.size(); pb++) {
               byphase[pb].hits += ss.phase_stats[pb].hits ;
               byphase[pb].moves += ss.phase_stats[pb].moves ;
               byphase[pb].expectedcost += ss.phase_stats[pb].expectedcost ;
            }
         }
      } ;
      vector<thread> threads ;
      for (int t = 0; t < nprob_threads; t++)
         threads.emplace_back(worker) ;
      for (auto &t : threads)
         t.join() ;
      r += total_mvs.load() ;
      cancellations += total_cancellations.load() ;
      blocks += total_blocks.load() ;
   } else {
      vector<int> xs ;
      vector<int> ys ;
      vector<ll> bits ;
      vector<int> next, prev ;
      vector<vector<planmask>> movemasks((size_t)phasewidth) ;
      for (ll i=0; i<(ll)phasewidth; i++)
         movemasks[(size_t)i].resize((size_t)cube.Mx) ;
      for (ll i=0; i<(ll)cube.My; i += (ll)phasewidth) {
         ll pw = (ll)phasewidth ;
         if (i + pw > (ll)cube.My)
            pw = (ll)cube.My - i ;
         for (ll wi=0; wi<pw; wi++)
            cube.getrow((ll)(i+wi)) ;
         int phasebase = 0 ;
         for (auto &f: reductionplan) {
            for (ll ii=0; ii<pw; ii++)
               for (ll j=0; j<(ll)cube.Mx; j++) {
                  movemasks[ii][j] = f.net->getplan(cube.rows[i+ii][j]) ;
                  if (fast && !writemoveinc)
                     f.net->doplan(movemasks[ii][j], cube.rows[i+ii][j]) ;
               }
            for (int pi=f.net->STAGE-1; pi>=0; pi--, phasebase++) {
               auto &alg = f.net->algs[pi] ;
               auto &wplan = f.plan[pi] ;
               wplan.plan((cubesizetype)pw) ;
               if (phasebase >= (int)byphase.size()) {
                  byphase.resize(phasebase+1) ;
                  byphase[phasebase].w1 = wplan.w1 ;
                  byphase[phasebase].w2 = wplan.w2 ;
                  byphase[phasebase].count1 = wplan.count1 ;
                  byphase[phasebase].count2 = wplan.count2 ;
               }
               byphase[phasebase].expectedcost += wplan.cost ;
               byphase[phasebase].ba = &(f.net->algs[pi]) ;
               if (showdetails) {
                  cout << "Alg " << phasebase ;
                  alg.show() ;
                  cout << endl ;
               }
               ll ii = 0 ;
               for (ll iii=0; iii<(ll)(wplan.count1 + wplan.count2); iii++) {
                  int w = iii < wplan.count1 ? wplan.w1 : wplan.w2 ;
                  bits.clear() ;
                  bits.resize((size_t)cube.Mx) ;
                  for (int wi=0; wi<w; wi++)
                     for (ll j=0; j<(ll)cube.Mx; j++)
                        if (movemasks[(size_t)(ii+wi)][(size_t)j].test(pi))
                           bits[(size_t)j] |= 1LL<<wi ;
                  prev.clear() ;
                  prev.resize(1LL<<w, -1) ;
                  next.resize((size_t)cube.Mx) ;
                  for (ll i=0; i<(ll)cube.Mx; i++) {
                     next[(size_t)i] = prev[bits[(size_t)i]] ;
                     prev[bits[(size_t)i]] = (int)i ;
                  }
                  for (ll bp=1; bp<(1LL<<w); bp++) {
                     if (prev[bp] < 0)
                        continue ;
                     ys.clear() ;
                     for (int wi=0; wi<w; wi++)
                        if ((bp >> wi) & 1)
                           ys.push_back(i+ii+wi) ;
                     xs.clear() ;
                     for (int i=prev[bp]; i>=0; i=next[i])
                        xs.push_back(i) ;
                     byphase[phasebase].hits += xs.size() * ys.size() ;
                     ll mvs = doblock(cube, alg, ys, xs) ;
                     if (writemoveinc && fast)
                        for (auto y: ys)
                           for (auto x: xs)
                              f.net->doplan(planmask{}.set(pi), cube.rows[y][x]) ;
                     byphase[phasebase].moves += mvs ;
                     r += mvs ;
                  }
                  ii += w ;
               }
            }
         }

         for (ll wi=0; wi<pw; wi++)
            cube.checkrelease((ll)(i+wi)) ;
         if (extrapolate)
            break ;
      }
   }
   if (dowings)
      for (ll k=0; k<cube.My; k++) {
         vector<int> w = cube.get_wings(k)[0] ;
         for (int t=0; t<24; t++)
            if (w[t] != t) {
               cerr << "wings NOT solved: orbit " << k << " slot " << t << endl ;
               exit(1) ;
            }
      }
   return r ;
}
ll dosolve(cubesizetype N, cubesizetype phasewidth) {
   prevmvs.clear() ;
   cubesizetype My = (N - 2) >> 1 ;
   cubesizetype Mx = (N - 2) - My ;
   blocks = 0 ;
   cancellations = 0 ;
   xcube cube(Mx, My, (readmoves || randmoves) ? 'I' : 'R') ;
   int fr[6] ; 
   for (int f=0; f<6; f++) { 
      fr[f] = cube.facerots[f] ; 
      cube.facerots[f] = 0 ; 
   }
   for (auto &row : cube.rows)
      for (auto &orb : row)
         for (int f=0; f<6; f++)
            rotate(orb.a + 4*f, orb.a + 4*f + (fr[f] & 3), orb.a + 4*f + 4) ;
   //cout << "=== wing orbits ===" << endl ;
   //cube.show_wings() ;
   //cout << "===================" << endl ;
#if WRITEFRAMES
   runningmovecount = 0 ;
   if (writemoveinc)
      writeframe(cube) ;
#endif
   ll r = solve(cube, phasewidth) ;
   sendblock() ; // flush out any remaining moves
   if (alwaysinmem)
      for (int i=0; i<cube.My; i++)
         cube.checkrow(i) ;
#if WRITEFRAMES
   if (writemoveinc) {
      writeframe(cube) ;
      cout << endl ;
   }
#endif
   return r ;
}
struct basicstat {
   basicstat(const char *n_) : name(n_), n(0), sx(0), sxx(0), hi(0), lo(0) {}
   void add(double x) {
      x *= extrapolatemul ;
      n++ ;
      if (n == 1)
         hi = lo = x ;
      sx += x ;
      sxx += x * x ;
      hi = max(hi, x) ;
      lo = min(lo, x) ;
   }
   void report() {
      double av = sx / n ;
      double stddev = sqrt(n * sxx - sx * sx) / n ;
      cout << name << " " << av << " stddev " << stddev << " range " << lo << " .. " << hi ;
   }
   void avreport() {
      double av = sx / n ;
      cout << name << " " << av ;
   }
   const char *name ;
   ll n ;
   double sx, sxx, hi, lo ;
} ;
basicstat mvstat("moves"), cancelstat("cancels"), blockstat("blocks") ;
void fullreport(ll at) {
   cout << "N " << N << " phasewidth " << phasewidth << " cnt " << at << " " ;
   mvstat.report() ; cout << " " ;
   cancelstat.avreport() ; cout << " " ;
   blockstat.avreport() ; cout << endl;
}
void addcounts(ll v) {
   mvstat.add(v)  ;
   cancelstat.add(cancellations) ;
   blockstat.add(blocks) ;
   if (solvecount > 1 && (mvstat.n & (mvstat.n - 1)) == 0 && mvstat.n != solvecount)
      fullreport(mvstat.n) ;
}
void usage() {
   cout << "Brobdicube (C) 2026 Tomas Rokicki." << endl <<
      endl <<
      "Usage:  brobdicube [options] size" << endl <<
      endl <<
      "Options:" << endl <<
      endl <<
      "-d      show some internal information on the strategy." << endl <<
      "-i      read scramble moves from stdin." << endl <<
      "-e      plan and estimate but do not solve." << endl <<
      "-r n    scramble with n random moves (printed as SCRAMBLE: ...)." << endl <<
      "-R      scramble with an auto-picked length (50*N) random moves." << endl <<
      "-M n    do this many solves and report aggregate stats." << endl <<
      "-o      write solution to stdout (on lines starting with a space)." << endl <<
      "-p n    set the phase chunking width." << endl <<
      "-s      actually (slowly) simulate moves; for testing." << endl <<
      "-x      extrapolate statistics from a single strip." << endl <<
      "-a x y  select first phase of length x and second phase of length y." << endl << 
      "-S      show stats for each phase separately." << endl <<
      "-g      use only greedy planning (skip optimal)." << endl <<
      "-show_wings print the edge permutations of the scrambled cube." << endl << 
#if WRITEFRAMES
      "-v n    write a png file of the cube state every n moves." << endl <<
#endif
      "-w      optimize plan for worst-case, not average, solves." << endl <<
      "-D dir  set the algs directory (default: ./algs3)." << endl <<
      "-t n    use n threads for probability cache precomputation (default: hardware_concurrency)." << endl <<
      "--describealgorithms" << endl <<
      "        for each loaded algorithm file, print (greedy and optimal separately):" << endl <<
      "        strategy, filename, depth (non-zero-prob algs), avg activations," << endl <<
      "        max activations, total cost (v/h/o/c 4-digit fields)," << endl <<
      "        and max single-alg cost.  No cube size required." << endl <<
      endl ;
   exit(0) ;
}
static string fmt_cost(ll cost) {
   ll v = cost / 1000000000000LL ;
   ll h = cost / 100000000LL % 10000LL ;
   ll o = cost / 10000LL % 10000LL ;
   ll c = cost % 10000LL ;
   char buf[32] ;
   snprintf(buf, sizeof(buf), "%04lld/%04lld/%04lld/%04lld", v, h, o, c) ;
   return buf ;
}
void describe_algorithms() {
   // Carry canonical key alongside each netgen so it can appear in output.
   vector<pair<string,netgen*>> all_ng ;
   for (auto& [ckey, gens] : algset)
      for (auto* ng : gens)
         all_ng.emplace_back(ckey, ng) ;
   sort(all_ng.begin(), all_ng.end(), [](const auto& a, const auto& b) {
      if (a.second->filename != b.second->filename)
         return a.second->filename < b.second->filename ;
      return (int)a.second->optimal < (int)b.second->optimal ; // greedy first
   }) ;
   int total = (int)all_ng.size() ;
   if (total == 0) return ;
   int nthreads = nprob_threads > 0 ? nprob_threads :
                  (int)thread::hardware_concurrency() ;
   if (nthreads < 1) nthreads = 1 ;
   if (nthreads > total) nthreads = total ;
   atomic<int> next_idx{0} ;
   auto worker = [&]() {
      for (;;) {
         int idx = next_idx.fetch_add(1, memory_order_relaxed) ;
         if (idx >= total) break ;
         auto& [canonical_key, ng] = all_ng[idx] ;
         // If max_activations_seen unknown (cache was loaded, not computed),
         // force a replan now to obtain it, then free the plan structures.
         if (ng->max_activations_seen < 0) {
            ng->plan() ;
            ng->plans.clear() ; ng->plans.shrink_to_fit() ;
         }
         double avg_act = 0 ;
         ll total_cost = 0, max_cost = 0 ;
         int depth = 0 ;
         for (int i = 0; i < (int)ng->algcosts.size(); i++) {
            if (ng->stageprob[i] > 0) {
               avg_act += ng->stageprob[i] ;
               depth++ ;
               total_cost += ng->algcosts[i] ;
               if (ng->algcosts[i] > max_cost)
                  max_cost = ng->algcosts[i] ;
            }
         }
         string mode = ng->optimal ? "optimal" : "greedy" ;
         lock_guard<mutex> lk(io_mutex) ;
         cout << mode
              << " key=" << ng->key
              << " file=" << ng->filename
              << " depth=" << depth
              << " avg_act=" << avg_act
              << " max_act=" << ng->max_activations_seen
              << " total_cost=" << fmt_cost(total_cost)
              << " max_alg_cost=" << fmt_cost(max_cost)
              << "\n" ;
      }
   } ;
   vector<thread> threads ;
   threads.reserve(nthreads) ;
   for (int i = 0; i < nthreads; i++)
      threads.emplace_back(worker) ;
   for (auto& t : threads)
      t.join() ;
}
void preload_probs() {
   vector<netgen*> all_ng ;
   for (auto& [key, gens] : algset)
      for (auto* ng : gens)
         all_ng.push_back(ng) ;
   int total = (int)all_ng.size() ;
   if (total == 0) return ;
   // Pre-scan: count netgens that need a cache computed (large enough to cache,
   // but cache file not yet present).  Only those will generate output.
   int need_compute = 0 ;
   for (auto* ng : all_ng) {
      if (ng->goal <= netgen::MINWRITE) continue ;
      FILE* f = ng->openprobcache("r") ;
      if (f) fclose(f) ; else need_compute++ ;
   }
   if (need_compute == 0) return ;  // all caches present — silent fast path
   int nthreads = nprob_threads > 0 ? nprob_threads :
                  (int)thread::hardware_concurrency() ;
   if (nthreads < 1) nthreads = 1 ;
   if (nthreads > total) nthreads = total ;
   cout << "Pre-computing probabilities for " << need_compute
        << " algorithm sets using " << nthreads << " threads" << endl ;
   atomic<int> next_idx{0} ;
   atomic<int> computed{0} ;
   auto worker = [&]() {
      for (;;) {
         int idx = next_idx.fetch_add(1, memory_order_relaxed) ;
         if (idx >= total) break ;
         bool was_computed = all_ng[idx]->ensureprobs() ;
         if (was_computed && all_ng[idx]->goal > netgen::MINWRITE) {
            int done = computed.fetch_add(1, memory_order_relaxed) + 1 ;
            lock_guard<mutex> lk(io_mutex) ;
            cout << "  probabilities: " << done << "/" << need_compute << " done" << endl ;
         }
      }
   } ;
   vector<thread> threads ;
   threads.reserve(nthreads) ;
   for (int i = 0; i < nthreads; i++)
      threads.emplace_back(worker) ;
   for (auto& t : threads)
      t.join() ;
   int bad_count = 0 ;
   for (auto ng : all_ng)
      if (ng->bad) bad_count++ ;
   if (bad_count > 0)
      cerr << "WARNING: " << bad_count << " inconsistent algorithm file(s) skipped." << endl ;
   cout << "All probability caches ready in " << duration() << endl ;
}
int main(int argc, char *argv[]) {
   duration() ;
   cout << "-" ;
   for (int i=0; i<argc; i++)
      cout << " " << argv[i] ;
   cout << endl ;
   while (argc > 1 && argv[1][0] == '-') {
      argc-- ;
      argv++ ;
      switch (argv[0][1]) {
case 'x':
         extrapolate = 1 ;
         break ;
case 'e':
         planonly = 1 ;
         break ;
case 's':
         fast = 0 ;
         break ;
case 'S':
         stats = 1 ;
         break ;
case 'g':
         greedyonly = 1 ;
         break ;
case 'd':
         showdetails = 1 ;
         break ;
case 'w':
         useavg = 0 ;
         break ;
case 'p':
         argc-- ;
         argv++ ;
         phasewidth = atocsz(argv[0]) ;
         break ;
case 'M':
         argc-- ;
         argv++ ;
         solvecount = atoll(argv[0]) ;
         break ;
#if WRITEFRAMES
case 'v':
         argc-- ;
         argv++ ;
         writemoveinc = atoll(argv[0]) ;
         break ;
#endif
case 'o':
         writemoves = 1 ;
         break ;
case 'i':
         readmoves = 1 ;
         break ;
case 'r':
         argc-- ;
         argv++ ;
         randmoves = atoll(argv[0]) ;
         break ;
case 'R':
         randmoves = -1 ;   // auto-pick length
         break ;
case 'P':
         argc-- ;
         argv++ ;
         plan_range_start = atocsz(argv[0]) ;
         argc-- ;
         argv++ ;
         plan_range_end = atocsz(argv[0]) ;
         break ;
case 'D':
         argc-- ;
         argv++ ;
         algsdir = argv[0] ;
         break ;
case 't':
         argc-- ;
         argv++ ;
         nprob_threads = atoi(argv[0]) ;
         break ;
default:
         if (strcmp(argv[0], "--describealgorithms") == 0) {
            describealgorithms = true ;
            break ;
         }
         error("! option not understood") ;
      }
   }
   if (!fast)
      nprob_threads = 1 ;
   expandrotations() ;
   build_wing_cycle_index() ;
   loadwingalgs() ;
   if (!fast && extrapolate)
      error("! cannot go slow and extrapolate at the same time") ;
   if (plan_range_start == 0 && !describealgorithms) {
      if (argc < 2)
         usage() ;
      if (argc > 2)
         error("! extra arguments given on command line") ;
      N = atocsz(argv[1]) ;
      if (randmoves < 0) // -1 means auto-pick a length of scramble
         randmoves = 50 * N ;
   }
   cout << setprecision(15) ;
   readalgs() ;
   preload_probs() ;
   if (describealgorithms) {
      describe_algorithms() ;
      return 0 ;
   }
   if (plan_range_start > 0) {
      const char *costname = useavg ? "average" : "worst-case" ;
      for (cubesizetype sz = plan_range_start; sz <= plan_range_end; sz *= 2) {
         cubesizetype My = (sz - 2) / 2 ;
         cubesizetype Mx = (sz - 2) - My ;
         cubesizetype pw = (phasewidth > 0) ? phasewidth : My ;
         bestcostcache.clear() ;
         bestcostplan.clear() ;
         reductionplan.clear() ;
         double cc = bestcost("000000", Mx, pw, useavg ? 'a' : 'w') ;
         if (pw != My)
            cc = cc * My / pw ;
         cout << "N=" << sz << " Estimated " << costname << " cost is " << cc << endl ;
         explain("000000") ;
      }
      exit(0) ;
   }
   cubesizetype My = (N - 2) / 2 ;
   cubesizetype Mx = (N - 2) - My ;
   if (phasewidth <= 0)
      phasewidth = My ;
   if (extrapolate)
      extrapolatemul = My / (double)phasewidth ;
   if (!fast || writemoveinc || randmoves || (readmoves && nprob_threads <= 1))
      alwaysinmem = 1 ;
   if ((readmoves || writemoves || writemoveinc) && solvecount > 1)
      error("! cannot write frames and read or write moves, or do multiple solves at the same time") ;
   double cc = bestcost("000000", Mx, phasewidth, useavg ? 'a' : 'w') ;
   if (randmoves)
      random_move_scramble(N, randmoves) ;
   else if (readmoves)
      readmoveset() ;
   else if (phasewidth != My)
      cc = cc * My / phasewidth ;
   cout << "Estimated " << (useavg ? "average" : "worst-case") << " cost is " << cc << endl ;
   reductionplan.clear() ;
   explain("000000") ;
   if (planonly)
      exit(0) ;
   ll mmc = 0 ;
   for (ll ii=1; ii<=solvecount; ii++) {
      ll t = dosolve(N, phasewidth) ;
      addcounts(t) ;
      mmc += t ;
   }
   if (solvecount == 1) {
      ll r = mmc ;
      if (extrapolate) {
         if (stats)
            for (int i=0; i<(int)byphase.size(); i++)
               byphase[i].showstats(i, My*phasewidth) ;
         lll lll_r = (lll)r * My / phasewidth ;
 cout << "r " << r << " My " << My << " lll_r " << lll_r << endl ;
         lll lll_cancellations = (lll)cancellations * My / phasewidth ;
         lll lll_blocks = (lll)blocks * My / phasewidth ;
         cout << "N " << N << " width " << phasewidth
              << " moves " << lll_r << " cancel " << lll_cancellations
              << " blocks " << lll_blocks << endl ;
      } else {
         if (stats)
            for (int i=0; i<(int)byphase.size(); i++)
               byphase[i].showstats(i, Mx*My) ;
         cout << "N " << N
              << " moves " << r << " cancel " << cancellations
              << " blocks " << blocks
              << " wingmoves " << wingmoves
              << " centermoves " << (r - wingmoves) << endl ;
      }
   } else {
      fullreport(solvecount) ;
   }
}
