/* Copyright 2015 Piotr Luszczek <luszczek@icl.utk.edu>
  
  Redistribution and use in source and binary forms, with or without
  modification,are permitted provided that the following conditions are met:
  
  1. Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
  
  2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
  
  3. Neither the name of the copyright holder nor the names of its contributors
  may be used to endorse or promote products derived from this software without
  specific prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* A C translation of the LAPACK fortran function DLAMCH(). Taken directly from
  http://icl.cs.utk.edu/lapack-forum/archives/lapack/msg01682.html
  with only trivial modifications (May 2017 Drew Schmidt). */


#ifndef __LINMOD2_MATLIB_DLAMCH__
#define __LINMOD2_MATLIB_DLAMCH__


#include <float.h>

#ifndef FLT_DIGITS
#define FLT_DIGITS 24
#endif
#ifndef DBL_DIGITS
#define DBL_DIGITS 53
#endif

static inline float slamch(const char *const restrict cmach)
{
  char ch = cmach[0];
  float sfmin, small, one = 1.0, zero = 0.0;
  
  if ('B' == ch || 'b' == ch) {
    return FLT_RADIX;
  } else if ('E' == ch || 'e' == ch) {
    return FLT_EPSILON;
  } else if ('L' == ch || 'l' == ch) {
    return FLT_MAX_EXP;
  } else if ('M' == ch || 'm' == ch) {
    return FLT_MIN_EXP;
  } else if ('N' == ch || 'n' == ch) {
    return FLT_DIGITS;
  } else if ('O' == ch || 'o' == ch) {
    return FLT_MAX;
  } else if ('P' == ch || 'p' == ch) {
    return FLT_EPSILON * FLT_RADIX;
  } else if ('R' == ch || 'r' == ch) {
    return FLT_ROUNDS < 2 ? one : zero;
  } else if ('S' == ch || 's' == ch) {
    /* Use SMALL plus a bit, to avoid the possibility of rounding causing overflow
      when computing  1/sfmin. */
    sfmin = FLT_MIN;
    small = one / FLT_MAX;
    if (small >= sfmin)
      sfmin = small * (one + FLT_EPSILON);
    return sfmin;
  } else if ('U' == ch || 'u' == ch) {
    return FLT_MIN;
  }

  return zero;
}



static inline double dlamch(const char *const restrict cmach)
{
  char ch = cmach[0];
  double sfmin, small, one = 1.0, zero = 0.0;

  if ('B' == ch || 'b' == ch) {
    return FLT_RADIX;
  } else if ('E' == ch || 'e' == ch) {
    return DBL_EPSILON;
  } else if ('L' == ch || 'l' == ch) {
    return DBL_MAX_EXP;
  } else if ('M' == ch || 'm' == ch) {
    return DBL_MIN_EXP;
  } else if ('N' == ch || 'n' == ch) {
    return DBL_DIGITS;
  } else if ('O' == ch || 'o' == ch) {
    return DBL_MAX;
  } else if ('P' == ch || 'p' == ch) {
    return DBL_EPSILON * FLT_RADIX;
  } else if ('R' == ch || 'r' == ch) {
    return FLT_ROUNDS < 2 ? one : zero;
  } else if ('S' == ch || 's' == ch) {
    /* Use SMALL plus a bit, to avoid the possibility of rounding causing overflow
      when computing  1/sfmin. */
    sfmin = DBL_MIN;
    small = one / DBL_MAX;
    if (small >= sfmin)
      sfmin = small * (one + DBL_EPSILON);
    return small;
  } else if ('U' == ch || 'u' == ch) {
    return DBL_MIN;
  }

  return zero;
}


#endif
