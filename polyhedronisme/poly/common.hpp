//
//  common.hpp
//  test_polygon
//
//  Created by asd on 03/09/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#pragma once

#include <stdlib.h>

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>
#include <string.h>

#include "simd.h"
#include "Timer.h"

using std::vector, std::string, std::map, std::to_string, std::pair,
    std::reverse, std::min, std::max, std::tuple, std::copy, std::get,
    std::lower_bound, std::inserter, std::unordered_set, std::atoi, std::atof;

using Vertex = Simd;  // simd_float3;
using Vertexes = vector<Vertex>;
using Face = vector<int>;
using Faces = vector<Face>;
using VertexesFloat = vector<vector<float>>;

class VertexIndex {
 public:
  int index;
  Vertex vertex;

  inline VertexIndex() : index(0), vertex() {}
  inline VertexIndex(int index, Vertex vertex) : index(index), vertex(vertex) {}
};

class Int4 {
 public:
  int i0, i1, i2, i3;

  inline Int4() {}
  inline Int4(int i0, int i1, int i2, int i3)
      : i0(i0), i1(i1), i2(i2), i3(i3) {}
  inline Int4(int i0, int i1, int i2) : i0(i0), i1(i1), i2(i2), i3(0) {}
  inline Int4(int i0, int i1) : i0(i0), i1(i1), i2(0), i3(0) {}
  inline Int4(int i0) : i0(i0), i1(0), i2(0), i3(0) {}

  inline const int &operator[](int i) const {
    switch (i) {
      default:
      case 0:
        return i0;
      case 1:
        return i1;
      case 2:
        return i2;
      case 3:
        return i3;
    }
  }  // accessor
  inline int &operator[](int i) {
    switch (i) {
      default:
      case 0:
        return i0;
      case 1:
        return i1;
      case 2:
        return i2;
      case 3:
        return i3;
    }
  }  // mutator

  inline bool operator<(const Int4 &o) const {
    if (i0 != o.i0) {
      return i0 < o.i0;
    } else {
      if (i1 != o.i1)
        return i1 < o.i1;
      else {
        if (i2 != o.i2)
          return i2 < o.i2;
        else
          return i3 < o.i3;
      }
    }
  }

  inline bool operator>(const Int4 &o) const {
    if (i0 != o.i0) {
      return i0 > o.i0;
    } else {
      if (i1 != o.i1)
        return i1 > o.i1;
      else {
        if (i2 != o.i2)
          return i2 > o.i2;
        else
          return i3 > o.i3;
      }
    }
  }

  inline bool operator==(const Int4 &o) const {
    return i0 == o.i0 && i1 == o.i1 && i2 == o.i2 && i3 == o.i3;
  }
  inline bool operator!=(const Int4 &o) const { return !(*this == o); }
};

class MapIndex {
 public:
  Int4 i0, i1, i2;

  inline MapIndex(Int4 i0, Int4 i1, Int4 i2) {
    this->i0 = i0;
    this->i1 = i1;
    this->i2 = i2;
  }
  inline MapIndex(Int4 i0, Int4 i1) { *this = MapIndex(i0, i1, 0); }
  inline MapIndex(Int4 i0) { *this = MapIndex(i0, 0, 0); }
  inline MapIndex() { *this = MapIndex(0, 0, 0); }

  static inline bool less(const MapIndex &a, const MapIndex &b) {
    return a < b;
  }

  inline bool operator<(const MapIndex &o) const {
    if (i0 < o.i0) return true;
    if (i0 > o.i0) return false;
    if (i1 < o.i1) return true;
    if (i1 > o.i1) return false;
    if (i2 < o.i2) return true;
    return false;  // i2>=o.i2
  }

  inline bool operator>(const MapIndex &o) const {
    if (i0 > o.i0) return true;
    if (i0 < o.i0) return false;
    if (i1 > o.i1) return true;
    if (i1 < o.i1) return false;
    if (i2 > o.i2) return true;
    return false;  // i2<=o.i2
  }
  inline const Int4 &operator[](int i) const {  // accessor
    switch (i) {
      case 0:
        return i0;
      case 1:
        return i1;
      case 2:
        return i2;
      default:
        return i0;
    }
  }
  inline Int4 &operator[](int i) {  // mutator
    switch (i) {
      case 0:
        return i0;
      case 1:
        return i1;
      case 2:
        return i2;
      default:
        return i0;
    }
  }
  inline bool operator==(const MapIndex &o) const {
    return memcmp(this, &o, sizeof(MapIndex)) == 0;
  }
  inline bool operator!=(const MapIndex &o) const {
    return memcmp(this, &o, sizeof(MapIndex)) != 0;
  }
};

class I4Vix {
 public:
  inline I4Vix() : index({}), vix() {}
  inline I4Vix(Int4 index, VertexIndex vix) : index(index), vix(vix) {}

  inline bool operator<(const I4Vix &o) { return index < o.index; }
  static inline bool less(const I4Vix &a, const Int4 &index) {
    return a.index < index;
  }
  inline bool operator==(const I4Vix &o) {
    return memcmp(this, &o, sizeof(*this)) == 0;
  }

  Int4 index;
  VertexIndex vix;
};

static inline Int4 to_int4(int v) { return {v + 1, 0, 0, 0}; }
static inline Int4 to_int4(int v1, int v2) { return {v1 + 1, v2 + 1, 0, 0}; }
static inline Int4 to_int4(int v1, int v2, int v3) {
  return {v1 + 1, v2 + 1, v3 + 1, 0};
}
static inline Int4 to_int4(int v1, int v2, int v3, int v4) {
  return {v1 + 1, v2 + 1, v3 + 1, v4 + 1};
}
static inline Int4 i4(int v1) { return to_int4(v1); }
static inline Int4 i4(int v1, int v2) { return to_int4(v1, v2); }
static inline Int4 i4(int v1, int v2, int v3) { return to_int4(v1, v2, v3); }
static inline Int4 i4(int v1, int v2, int v3, int v4) {
  return to_int4(v1, v2, v3, v4);
}

static inline Int4 i4_min(int v1, int v2) {
  return v1 < v2 ? i4(v1, v2) : i4(v2, v1);
}
static inline Int4 i4_min(int i, int v1, int v2) {
  return v1 < v2 ? i4(i, v1, v2) : i4(i, v2, v1);
}

static string str(size_t i) { return to_string(i); }
static string str(int i) { return to_string(i); }
static inline string str(string s, size_t i) { return s + str(i); }
static inline string midName(int v1, int v2) {
  return v1 < v2 ? str(v1) + "_" + str(v2) : str(v2) + "_" + str(v1);
}

static Vertex calcCentroid(Vertexes vertices) {
  // running sum of vertex coords
  Vertex centroidV = 0;
  for (auto &v : vertices) centroidV += v;
  return centroidV / vertices.size();
}

static inline Vertex midpoint(Vertex vec1, Vertex vec2) {
  return (vec1 + vec2) / 2.;
}

static inline Vertex tween(Vertex &vec1, Vertex &vec2, float t) {
  return (vec1 * (1.f - t)) + (vec2 * t);
}
static inline Vertex oneThird(Vertex &vec1, Vertex &vec2) {
  return tween(vec1, vec2, 1 / 3.f);
}
static inline int intersect(Face &set1, Face &set2, Face &set3) {
  for (auto s1 : set1)
    for (auto s2 : set2)
      if (s1 == s2)
        for (auto s3 : set3)
          if (s1 == s3) return s1;
  return -1;  // empty intersection
}
