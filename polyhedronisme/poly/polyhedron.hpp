//
//  polyhedron.hpp
//  test_polyhedron
//
//  Created by asd on 03/09/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#pragma once

#include "Thread.h"
#include "color.hpp"
#include "common.hpp"

class Polyhedron {
 public:
  Polyhedron() {}

  Polyhedron(const string name, const Vertexes vertexes, const Faces faces)
      : name(name), vertexes(vertexes), faces(faces) {}

  Polyhedron(const string name, int c, const VertexesFloat vertexes,
             const vector<vector<int>> faces)
      : name(name), faces(faces) {
    for (auto v : vertexes) this->vertexes.push_back(Vertex{v[0], v[1], v[2]});
  }

  Polyhedron recalc() {
    calc_normals();
    calc_areas();
    calc_centers();
    calc_colors();
    return *this;
  }

  inline size_t n_vertexes() { return vertexes.size(); }
  inline size_t n_faces() { return faces.size(); }

  void scale_vertexes() {
    float min = __FLT_MAX__, max = -__FLT_MAX__;
    for (auto &v : vertexes) {
      max = std::max(max, simd_reduce_max(v));
      min = std::min(min, simd_reduce_min(v));
    }
    float diff = abs(max - min);
    if (diff != 0.f)
      for (auto &v : vertexes) v /= diff;
  }

  void calc_normals() {  // per face
    normals = Vertexes(n_faces());
    Thread(n_faces()).run([this](int f) {
      normals[f] = calc_normal(vertexes[faces[f][0]], vertexes[faces[f][1]],
                               vertexes[faces[f][2]]);
    });
  }

  void calc_normals_st() {  // per face
    normals = Vertexes(n_faces());
    for (size_t f = 0; f < n_faces(); f++)
      normals[f] = calc_normal(vertexes[faces[f][0]], vertexes[faces[f][1]],
                               vertexes[faces[f][2]]);
  }

  int count_points() {  // count # of vertex used in all faces
    int np = 0;
    for (auto &face : faces) np += face.size();
    return np;
  }

  // calculate average normal vector for array of vertices
  Vertexes avg_normals() {
    auto normals = Vertexes(n_faces());

    Thread(n_faces()).run([this, &normals](int i) {
      size_t face_len = faces[i].size();
      Vertex normalV = 0;
      auto v1 = vertexes[face_len - 2], v2 = vertexes[face_len - 1];

      for (auto ic : faces[i]) {  // running sum of normal vectors
        auto v3 = vertexes[ic];
        normalV += normal(v1, v2, v3);
        v1 = v2;
        v2 = v3;  // shift over one
      }
      normals[i] = unit(normalV);
    });
    return normals;
  }

  Vertexes avg_normals_st() {
    auto normals = Vertexes(n_faces());

    for (size_t i = 0; i < n_faces(); i++) {
      size_t face_len = faces[i].size();
      Vertex normalV = 0;
      auto v1 = vertexes[face_len - 2], v2 = vertexes[face_len - 1];

      for (auto ic : faces[i]) {  // running sum of normal vectors
        auto v3 = vertexes[ic];
        normalV += normal(v1, v2, v3);
        v1 = v2;
        v2 = v3;  // shift over one
      }
      normals[i] = unit(normalV);
    }
    return normals;
  }

  void calc_centers() {  // per face
    centers = Vertexes(n_faces());
    Thread(n_faces()).run([this](int f) {
      Vertex fcenter = 0;
      Face face = faces[f];
      // average vertex coords
      for (size_t ic = 0; ic < face.size(); ic++) fcenter += vertexes[face[ic]];
      centers[f] =
          fcenter / face.size();  //  return face - ordered array  of  centroids
    });
  }

  void calc_centers_st() {  // per face
    centers = Vertexes(n_faces());
    for (size_t f = 0; f < n_faces(); f++) {
      Vertex fcenter = 0;
      Face face = faces[f];
      // average vertex coords
      for (size_t ic = 0; ic < face.size(); ic++) fcenter += vertexes[face[ic]];
      centers[f] =
          fcenter / face.size();  //  return face - ordered array  of  centroids
    }
  }

  void calc_areas() {  // per face, requires normals
    areas = vector<float>(n_faces());
    Thread(n_faces()).run([this](int f) {
      auto &face = faces[f];
      Vertex vsum = 0;
      auto fl = face.size();
      Vertex v1 = vertexes[face[fl - 2]], v2 = vertexes[face[fl - 1]];

      for (int ic = 0; ic < fl; ic++) {
        vsum += cross(v1, v2);
        v1 = v2;
        v2 = vertexes[face[ic]];
      }
      areas[f] = abs(dot(normals[f], vsum)) / 2;
    });
  }

  void calc_areas_st() {  // per face
    areas = vector<float>(n_faces());
    for (size_t f = 0; f < n_faces(); f++) {
      auto &face = faces[f];
      Vertex vsum = 0;
      auto fl = face.size();
      Vertex v1 = vertexes[face[fl - 2]], v2 = vertexes[face[fl - 1]];

      for (size_t ic = 0; ic < fl; ic++) {
        vsum += cross(v1, v2);
        v1 = v2;
        v2 = vertexes[face[ic]];
      }
      areas[f] = abs(dot(normals[f], vsum)) / 2;
    }
  }

  void calc_colors() {  // per areas
    calc_areas();       // areas required

    auto pallette = Color::random_pallete();
    map<int, Vertex> color_dict;  // color dict<sigfigs, pallette>

    for (auto a : areas)
      color_dict.try_emplace(sigfigs(a),
                             pallette[color_dict.size() % pallette.size()]);
    colors.clear();
    for (auto a : areas) colors.push_back(color_dict[sigfigs(a)]);
  }

  Vertex centroid(Face &face) {
    Vertex centroid = 0;  // calc centroid of face
    for (auto ic : face) centroid += vertexes[ic];
    return centroid /= face.size();
  }

  // sets
  //  void set_faces(Faces &faces) {
  //    this->faces = faces;
  //    this->n_faces = faces.size();
  //  }

  // gets
  string get_name() { return name; }
  Vertexes get_vertexes() { return vertexes; }
  Faces get_faces() { return faces; }
  Vertexes get_normals() {
    if (normals.empty()) calc_normals();
    return normals;
  }
  vector<float> get_areas() {
    if (areas.empty()) calc_areas();
    return areas;
  }
  Vertexes get_centers() {
    if (centers.empty()) calc_centers();
    return centers;
  }
  Vertexes get_colors() {
    get_areas();  // required
    if (colors.empty()) calc_colors();
    return colors;
  }

  Vertex get_color(int face) {  // get face color according to face area
    return colors[face];
  }

  Vertex get_normal(int face) { return normals[face]; }

  void set_colors(Vertexes &colors) { this->colors = colors; }

  void set_normals(Vertexes &normals) { this->normals = normals; }

  void set_vertexes(Vertexes &vertexes) {
    this->vertexes = vertexes;
    
  }

  void replace(Vertexes vertexes, Vertexes normals, Vertexes colors) {
    this->vertexes = vertexes;
    
    this->normals = normals;
    this->colors = colors;
  }

  void new_colors() {
    colors.clear();
    calc_colors();
  }

  Polyhedron normalize()  // remove unused vertexes
  {
    vector<int> oldNew(maxFaceIndex() + 1);
    Face face;
    int nv = 0, nvdx = 0;
    Vertexes usedVtx;

    fill(oldNew.begin(), oldNew.end(), -1);  // oldNew = -1

    for (auto face : faces) {
      for (auto ix : face) {
        if (oldNew[ix] == -1) {
          oldNew[ix] = nvdx++;
          usedVtx.push_back(vertexes[ix]);
        }
      }
    }

    for (auto ix = 0; ix < faces.size(); ix++)  // assign faces
      for (int i = 0; i < faces[ix].size(); i++)
        faces[ix][i] = oldNew[faces[ix][i]];

    set_vertexes(usedVtx);  // assign vertexes & faces

    centers.clear();
    normals.clear();
    colors.clear();
    areas.clear();

    return *this;
  }
  // print
  void print_stat() {
    recalc();

    printf("name    : %s\n", name.c_str());
    printf("vertices: %ld\n", vertexes.size());
    printf("faces   : %ld\n", faces.size());

    for (int i = 0; i < faces.size(); i++) {
      auto face = faces[i];  // traverse vertexes
      for (auto ixc : face)
        printf("(%f, %f, %f) ", vertexes[ixc].v[0], vertexes[ixc].v[1],
               vertexes[ixc].v[2]);
    }
    puts("traverse ok");
  }

  string formatString(const char *format, ...) {
    char buffer[1024];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);
    return std::string(buffer);
  }

  string sprint() {
    recalc();

    string res;
    res += formatString("name:%s\n", name.c_str());

    res += formatString("vertices %ld\n", vertexes.size());
    int iv = 0;
    for (auto v : vertexes)
      res +=
          formatString("%d:(%.2f, %.2f, %.2f), ", iv++, v.v[0], v.v[1], v.v[2]);

    res += formatString("\nfaces %ld\n", faces.size());
    int iface = 0;
    for (int i = 0; i < faces.size(); i++) {
      auto face = faces[i];
      res += formatString("%d:(", iface++);
      for (auto ic : face) res += formatString("%d ", ic);
      res += formatString("), ");
    }

    res += formatString("\nnormals\n");
    for (size_t i = 0; i < faces.size(); i++)
      res += formatString("%ld:(%.2f, %.2f, %.2f), ", i, normals[i].v[0],
                          normals[i].v[1], normals[i].v[2]);

    res += formatString("\nareas\n");
    for (size_t i = 0; i < faces.size(); i++)
      res += formatString("%ld:(%.2f), ", i, areas[i]);

    res += formatString("\ncolors\n");
    for (size_t i = 0; i < colors.size(); i++)
      res += formatString("%ld:(%.2f, %.2f, %.2f, %s), ", i, colors[i].v[0],
                          colors[i].v[1], colors[i].v[2],
                          Color::rgb2hex(colors[i]).c_str());

    res += formatString("\ncenters\n");
    for (size_t i = 0; i < faces.size(); i++)
      res += formatString("%ld:(%.2f, %.2f, %.2f), ", i, centers[i].v[0],
                          centers[i].v[1], centers[i].v[2]);
    res += formatString("\n");

    return res;
  }

  void print() {
    recalc();

    printf("name:%s\n", name.c_str());

    printf("vertices %ld\n", vertexes.size());
    int iv = 0;
    for (auto v : vertexes)
      printf("%d:(%.2f, %.2f, %.2f), ", iv++, v.x(), v.y(), v.z());

    printf("\nfaces %ld\n", faces.size());
    int iface = 0;
    for (int i = 0; i < faces.size(); i++) {
      auto face = faces[i];
      printf("%d:(", iface++);
      for (auto ic : face) printf("%d ", ic);
      printf("), ");
    }

    puts("\nnormals");
    for (size_t i = 0; i < faces.size(); i++)
      printf("%ld:(%.2f, %.2f, %.2f), ", i, normals[i].x(), normals[i].y(),
             normals[i].z());

    puts("\nareas");
    for (size_t i = 0; i < faces.size(); i++)
      printf("%ld:(%.2f), ", i, areas[i]);

    puts("\ncolors");
    for (size_t i = 0; i < colors.size(); i++)
      printf("%ld:(%.2f, %.2f, %.2f, %s), ", i, colors[i].x(), colors[i].y(),
             colors[i].z(), Color::rgb2hex(colors[i]).c_str());

    puts("\ncenters");
    for (size_t i = 0; i < faces.size(); i++)
      printf("%ld:(%.2f, %.2f, %.2f), ", i, centers[i].x(), centers[i].y(),
             centers[i].z());
    puts("");
  }

 public:
  string name = "";
  Vertexes vertexes;
  Faces faces;
  // size_t n_vertex = 0, n_faces = 0;

 public:
  Vertexes normals, colors, centers;
  vector<float> areas;

 private:
  inline Vertex calc_normal(Vertex v0, Vertex v1, Vertex v2) {
    return unit(cross(v1 - v0, v2 - v1));
  }
  inline Vertex normal(Vertex v0, Vertex v1, Vertex v2) {
    return cross(v1 - v0, v2 - v1);
  }

  // inline Vertex unit(const Vertex v) { return normalize(v); }

  int maxFaceIndex() {
    int res = -1;
    for (auto face : faces)
      for (auto ix : face) res = max(res, ix);
    return res;
  }

  int nVertexes() {
    int res = 0;
    for (auto face : faces) res += face.size();
    return res;
  }

  int sigfigs(float f,
              int nsigs = 2) {  // returns w. nsigs digits ignoring magnitude
    if (f == 0.f) return 0;
    float mantissa = f / powf(10, floor(log10f(f)));
    return int(roundf(mantissa * powf(10, (nsigs - 1))));
  }
  int sigfigs_fast(
      float f) {  // returns string w. nsigs digits ignoring magnitude
    return int(f * 100);
  }
};

