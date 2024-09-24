#pragma once

#include <math.h>
using std::max;
using std::min;

class Simd {
  typedef float simd_f3 __attribute__((ext_vector_type(3)));
  typedef int v4si __attribute__((vector_size(16)));

 public:
  simd_f3 v;

 public:
  inline Simd() { v = {0}; }
  inline Simd(simd_f3 v) : v(v) {}
  // inline Simd(Simd &v) : v(v.v) {}
  // inline Simd(float f) : v(f) {}
  inline Simd(float r, float g, float b) {
    v[0] = r;
    v[1] = g;
    v[2] = b;
  }
  inline float x() { return v[0]; }
  inline float y() { return v[1]; }
  inline float z() { return v[2]; }

  inline Simd operator+(Simd vo) { return v + vo.v; }
  inline Simd operator-(Simd vo) { return v - vo.v; }
  inline Simd operator*(Simd vo) { return v * vo.v; }
  inline Simd operator*(float f) { return f * v; }

  inline Simd operator/(Simd vo) { return v / vo.v; }
  inline Simd operator-() { return -v; }

  inline Simd operator+=(Simd vo) {
    v += vo.v;
    return *this;
  }
  inline Simd operator-=(Simd vo) {
    v -= vo.v;
    return *this;
  }
  inline Simd operator*=(Simd vo) {
    v *= vo.v;
    return *this;
  }
  inline Simd operator/=(Simd vo) {
    v /= vo.v;
    return *this;
  }

  inline float len_sq() { return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; }
  inline float len() { return sqrtf(len_sq()); }
  inline float sum() { return v[0] + v[1] + v[2]; }
  inline Simd unit() {
    float l = len();
    return (l != 0.0f) ? v / len() : v;
  }
};

inline Simd mid(Simd v1, Simd v2) { return (v1 + v2) / 2.0f; }
inline float len_sq(Simd v) { return v.len_sq(); }
inline float len(Simd v) { return v.len(); }
inline float dot(Simd v1, Simd v2) { return (v1 * v2).sum(); }
inline Simd unit(Simd v) { return v.unit(); }
static inline Simd cross(Simd v1, Simd v2) {
  return (v1.v.zxy * v2.v - v1.v * v2.v.zxy).zxy;
  // return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y *
  // b.x};
}
static inline Simd normalize(Simd v) { return v / len(v); }
static inline float simd_reduce_max(Simd v) {
  return max(v.v[0], max(v.v[1], v.v[2]));
}
static inline float simd_reduce_min(Simd v) {
  return min(v.v[0], max(v.v[1], v.v[2]));
}
static inline Simd operator*(float f, Simd vo) { return vo * f; }
