//
//  poly_operations.hpp
//  test_polygon
//
//  Created by asd on 05/09/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#pragma once
#pragma GCC diagnostic ignored "-Wmultichar"

#include "common.hpp"

#include "Thread.h"
#include "fastflags.h"
#include "polyhedron.hpp"

//===================================================================================================
// Polyhedron Operators
//===================================================================================================
// for each vertex of new polyhedron:
//     call newV(Vname, xyz) with a symbolic name and coordinates
// for each flag of new polyhedron:
//     call newFlag(Fname, Vname1, Vname2) with a symbolic name for the new face
//     and the symbolic name for two vertices forming an oriented edge
// ORIENTATION -must- be dealt with properly to make a manifold (correct) mesh.
// Specifically, no edge v1->v2 can ever be crossed in the -same direction- by
// two different faces
//
// call topoly() to assemble flags into polyhedron structure by following the
// orbits of the vertex mapping stored in the flagset for each new face
//
// set name as appropriate

class PolyOperations {
 public:
  // Kis(N)
  // ------------------------------------------------------
  // Kis (abbreviated from triakis) transforms an N-sided face into an N-pyramid
  // rooted at the same base vertices. only kis n-sided faces, but n==0 means
  // kis all.
  //
  static Polyhedron kisN(Polyhedron &poly, int n = 0, float apexdist = 0.1f) {
    int nth = Thread::getnthreads();
    vector<Flag> flags(nth);

    auto normals = poly.get_normals();
    auto centers = poly.get_centers();

    bool foundAny = false;

    // create face map
    Thread(poly.n_faces())
        .run([&flags, &poly, &foundAny, n, apexdist, &centers, &normals](
                 int t, int nface) {
          Flag &flag = flags[t];
          auto &face = poly.faces[nface];
          auto fname = i4('f', nface);

          int v1 = face.back();

          for (auto v2 : face) {
            auto iv2 = i4(v2);

            flag.add_vertex(iv2, poly.vertexes[v2]);  // poly.vtx

            if (face.size() == n || n == 0) {
              foundAny = true;

              flag.add_vertex(fname,
                              centers[nface] + (normals[nface] *
                                                apexdist));  // raised center
              flag.add_face({i4(v1), iv2, fname});
            } else {
              flag.add_face(nface, i4(v1), i4(v2));
            }

            v1 = v2;  // current becomes previous
          }
        });

    // combine flags->flag
    Flag flag(flags);

    if (!flag.valid)
      return poly;
    else
      return Polyhedron("k" + (n ? str(n) : "") + poly.name, flag.vertexes,
                        flag.faces)
          .normalize();
  }

  // Ambo
  // ------------------------------------------------------
  // The best way to think of the ambo operator is as a topological "tween"
  // between a polyhedron and its dual polyhedron.  Thus the ambo of a dual
  // polyhedron is the same as the ambo of the original. Also called
  // "Rectify".
  //

  static Polyhedron ambo(Polyhedron &poly) {
    int nth = Thread::getnthreads();
    vector<Flag> flags(nth);

    Thread(poly.n_faces()).run([&flags, &poly](int t, int nface) {
      Flag &flag = flags[t];
      auto &face = poly.faces[nface];
      auto flen = face.size();

      auto v1 = face[flen - 2],
           v2 = face[flen - 1];  //  [v1, v2,f.slice(-2);

      vector<Int4> f_orig;
      for (auto v3 : face) {
        auto m12 = i4_min(v1, v2), m23 = i4_min(v2, v3);

        if (v1 <
            v2)  // vertices are the midpoints of all edges of original poly
          flag.add_vertex(m12, midpoint(poly.vertexes[v1], poly.vertexes[v2]));

        // two new flags:
        // One whose face corresponds to the original f:
        f_orig.push_back(m12);

        // Another flag whose face  corresponds to (the truncated) v2:
        flag.add_face(i4('dual', v2), m23, m12);

        // shift over one
        v1 = v2;
        v2 = v3;
      }
      flag.add_face(f_orig);
    });

    Flag flag(flags);
    if (!flag.valid)
      return poly;
    else

      return Polyhedron("a" + poly.name, flag.vertexes, flag.faces).normalize();
  }

  // Gyro
  // -------------------------------------------------------
  // This is the dual operator to "snub", i.e dual*Gyro = Snub.  It is a bit
  // easier to implement this way.
  //
  // Snub creates at each vertex a new face, expands and twists it, and adds
  // two new triangles to replace each edge.

  static Polyhedron gyro(Polyhedron &poly) {
    Vertexes centers =
        poly.get_centers();  // new vertices in center of each face

    int nth = Thread::getnthreads();  // one flag per thread
    vector<Flag> flags(nth);

    Flag flag(poly.vertexes);

    Thread(poly.n_faces()).run([&flags, &poly, &centers](int t, int i) {
      Flag &flag = flags[t];
      auto &f = poly.faces[i];
      auto flen = f.size();

      auto v1 = f[flen - 2], v2 = f[flen - 1];  //  [v1, v2,f.slice(-2);

      flag.add_vertex(i4('cntr', i), centers[i]);

      for (size_t j = 0; j < flen; j++) {
        auto sv1 = str(v1), sv2 = str(v2), si = str(i);
        auto v = f[j];
        auto v3 = v;
        auto sv3 = str(v3);

        flag.add_vertex(
            i4(v1, v2),
            oneThird(poly.vertexes[v1], poly.vertexes[v2]));  // new v in face

        // 5 new faces
        flag.add_face(
            {i4('cntr', i), i4(v1, v2), i4(v2, v1), i4(v2), i4(v2, v3)});

        // shift over one
        v1 = v2;
        v2 = v3;
      }
    });

    flag.combine(flags);

    if (!flag.valid)
      return poly;
    else
      return Polyhedron("g" + poly.name, flag.vertexes, flag.faces).normalize();
  }

  // Propellor
  // ------------------------------------------------------
  // builds a new 'skew face' by making new points along edges, 1/3rd the
  // distance from v1->v2, then connecting these into a new inset face.  This
  // breaks rotational symmetry about the faces, whirling them into gyres

  static Polyhedron propellor(Polyhedron &poly) {
    Flag flag(poly.vertexes);
    vector<Flag> flags(Thread::getnthreads());  // one flag per thread

    Thread(poly.n_faces()).run([&flags, &poly](int t, int i) {
      Flag &flag = flags[t];
      auto &f = poly.faces[i];
      auto flen = f.size();
      auto v1 = f[flen - 2], v2 = f[flen - 1];  //  [v1, v2,f.slice(-2);

      for (auto v3 : f) {
        flag.add_vertex(
            i4(v1, v2),
            oneThird(poly.vertexes[v1],
                     poly.vertexes[v2]));  // new v in face, 1/3rd along edge

        flag.add_face(i4(i), i4(v1, v2), i4(v2, v3));  // five new flags
        flag.add_face({i4(v1, v2), i4(v2, v1), i4(v2), i4(v2, v3)});

        // shift over one
        v1 = v2;
        v2 = v3;
      }
    });

    flag.combine(flags);

    if (!flag.valid)
      return poly;
    else
      return Polyhedron("p" + poly.name, flag.vertexes, flag.faces).normalize();
  }

  // Reflection
  // ------------------------------------------------------
  // geometric reflection through origin

  static Polyhedron reflect(Polyhedron &poly) {
    // reflect each point through origin
    Thread(poly.n_vertexes()).run([&poly](int i) {
      poly.vertexes[i] = -poly.vertexes[i];
    });

    // repair clockwise-ness of faces
    Thread(poly.n_faces()).run([&poly](int i) {
      reverse(poly.faces[i].begin(), poly.faces[i].end());
    });

    poly.name = "r" + poly.name;
    return poly;
  }

  // Dual
  // ---------------------------------------------------------
  // The dual of a polyhedron is another mesh wherein:
  // - every face in the original becomes a vertex in the dual
  // - every vertex in the original becomes a face in the dual
  //
  // So N_faces, N_vertices = N_dualfaces, N_dualvertices
  //
  // The new vertex coordinates are convenient to set to the original face
  // centroids.
  //
  static Polyhedron dual(Polyhedron &poly) {
    auto face_map = Flag::gen_face_map(poly);
    auto centers = poly.get_centers();
    vector<Flag> flags(Thread::getnthreads());  // one flag per thread

    Thread(poly.n_faces())
        .run([&flags, &centers, &face_map, &poly](int t, int i) {
          Flag &flag = flags[t];
          auto &f = poly.faces[i];
          auto v1 = f.back();  // previous vertex
          flag.add_vertex(i4(i), centers[i]);
          for (auto v2 : f) {
            flag.add_face(i4(v1), Flag::Int4int::find(face_map, i4(v2, v1)),
                          i4(i));
            v1 = v2;  // current becomes previous
          }
        });

    auto &pn = poly.name;
    auto name = (pn[0] != 'd') ? "d" + pn : pn.substr(1, string::npos);

    Flag flag(flags);

    if (!flag.valid)
      return poly;
    else
      return Polyhedron(name, flag.vertexes, flag.faces).normalize();
  }

  // Chamfer
  // ----------------------------------------------------
  // A truncation along a polyhedron's edges.
  // Chamfering or edge-truncation is similar to expansion, moving faces apart
  // and outward, but also maintains the original vertices. Adds a new
  // hexagonal face in place of each original edge. A polyhedron with e edges
  // will have a chamfered form containing 2e new vertices, 3e new edges, and
  // e new hexagonal faces. -- Wikipedia See also
  // http://dmccooey.com/polyhedra/Chamfer.html
  //
  // The dist parameter could control how deeply to chamfer.
  // But I'm not sure about implementing that yet.
  //
  // Q: what is the dual operation of chamfering? I.e.
  // if cX = dxdX, and xX = dcdX, what operation is x?

  // We could "almost" do this in terms of already-implemented operations:
  // cC = t4daC = t4jC, cO = t3daO, cD = t5daD, cI = t3daI
  // But it doesn't work for cases like T.

  static Polyhedron chamfer(Polyhedron &poly, float dist = 0.05) {
    int nth = Thread::getnthreads();
    vector<Flag> flags(nth);
    auto normals = poly.get_normals();

    // For each face f in the original poly
    Thread(poly.n_faces()).run([&poly, &flags, dist, &normals](int t, int i) {
      auto &flag = flags[t];
      auto &f = poly.faces[i];
      auto v1 = f.back();
      auto v1new = i4(i, v1);

      for (auto &v2 : f) {
        // TODO: figure out what distances will give us a planar hex face.
        // Move each old vertex further from the origin.
        flag.add_vertex(i4(v2), (1.0f + dist) * poly.vertexes[v2]);
        // Add a new vertex, moved parallel to normal.
        auto v2new = i4(i, v2);

        flag.add_vertex(v2new, poly.vertexes[v2] + (dist * 1.5f * normals[i]));

        // Four new flags:
        // One whose face corresponds to the original face:
        flag.add_face(i4('orig', i), v1new, v2new);

        // And three for the edges of the new hexagon:
        auto facename = v1 < v2 ? i4('hex', v1, v2) : i4('hex', v2, v1);
        flag.add_face(facename, i4(v2), v2new);
        flag.add_face(facename, v2new, v1new);
        flag.add_face(facename, v1new, i4(v1));

        v1 = v2;
        v1new = v2new;
      }
    });

    Flag flag(flags);

    if (!flag.valid)
      return poly;
    else
      return Polyhedron("c" + poly.name, flag.vertexes, flag.faces).normalize();
  }

  // Whirl
  // -------------------------------------------------------
  // Gyro followed by truncation of vertices centered on original faces.
  // This create 2 new hexagons for every original edge.
  // (https://en.wikipedia.org/wiki/Conway_polyhedron_notation#Operations_on_polyhedra)
  //
  // Possible extension: take a parameter n that means only whirl n-sided
  // faces. If we do that, the flags marked #* below will need to have their
  // other sides filled in one way or another, depending on whether the
  // adjacent face is whirled or not.

  static Polyhedron whirl(Polyhedron &poly, int n = 0) {
    (void)n;

    int nth = Thread::getnthreads();
    vector<Flag> flags(nth);
    Flag flag(poly.vertexes);

    // new vertices around center of each face
    auto centers = poly.get_centers();

    Thread(poly.n_faces()).run([&poly, &flags, &centers](int t, int i) {
      auto &flag = flags[t];
      auto &f = poly.faces[i];
      auto flen = f.size();
      auto v1 = f[flen - 2], v2 = f[flen - 1];  //  [v1, v2,f.slice(-2);

      for (size_t j = 0; j < flen; j++) {
        auto v = f[j];
        auto v3 = v;

        // New vertex along edge
        auto v1_2 = oneThird(poly.vertexes[v1], poly.vertexes[v2]);
        flag.add_vertex(i4(v1, v2), v1_2);
        // New vertices near center of face

        auto cv1name = i4('cntr', i, v1);
        auto cv2name = i4('cntr', i, v2);

        flag.add_vertex(cv1name, unit(oneThird(centers[i], v1_2)));

        //        auto fname = i4(i, 'f', v1);
        // New hexagon for each original edge
        flag.add_face(
            {cv1name, i4(v1, v2), i4(v2, v1), i4(v2), i4(v2, v3), cv2name});

        // New face in center of each old face
        flag.add_face(i4('c', i), cv1name, cv2name);

        v1 = v2;  // shift over one
        v2 = v3;
      }
    });

    flag.combine(flags);

    if (!flag.valid)
      return poly;
    else
      return Polyhedron("w" + poly.name, flag.vertexes, flag.faces).normalize();
  }

  // Quinto
  // -------------------------------------------------------
  // This creates a pentagon for every point in the original face, as well as
  // one new inset face.
  static Polyhedron quinto(Polyhedron &poly) {
    vector<Flag> flags(Thread::getnthreads());

    auto centers = poly.get_centers();

    Thread(poly.n_faces()).run([&flags, &poly, &centers](int t, int nface) {
      Flag &flag = flags[t];

      // For each face f in the original poly

      auto &f = poly.faces[nface];
      auto flen = f.size();
      auto &centroid = centers[nface];

      // walk over face vertex-triplets
      auto v1 = f[flen - 2], v2 = f[flen - 1];  //  [v1, v2,f.slice(-2);

      vector<Int4> vi4;
      for (auto v3 : f) {
        auto t12 = i4_min(v1, v2), ti12 = i4_min(nface, v1, v2),
             t23 = i4_min(v2, v3), ti23 = i4_min(nface, v2, v3), iv2 = i4(v2);

        // for each face-corner, we make two new points:
        Vertex midpt = midpoint(poly.vertexes[v1], poly.vertexes[v2]),
               innerpt = midpoint(midpt, centroid);

        flag.add_vertex(t12, midpt);
        flag.add_vertex(ti12, innerpt);

        // and add the old corner-vertex
        flag.add_vertex(iv2, poly.vertexes[v2]);

        // pentagon for each vertex in original face

        flag.add_face({ti12, t12, iv2, t23, ti23});

        // inner rotated face of same vertex-number as original
        vi4.push_back(ti12);

        // shift over one
        v1 = v2;
        v2 = v3;
      }
      flag.add_face(vi4);
    });

    //    Flag flag(flags);
    Flag flag;
    flag.combine(flags);

    if (!flag.valid)
      return poly;
    else
      return Polyhedron("q" + poly.name, flag.vertexes, flag.faces).normalize();
  }

  // inset / extrude / "Loft" operator
  // ------------------------------------------------------
  static Polyhedron insetN(Polyhedron &poly, int n = 0, float inset_dist = 0.3f,
                           float popout_dist = -0.1f) {
    Flag flag(poly.vertexes);
    vector<Flag> flags(Thread::getnthreads());

    auto normals = poly.get_normals();
    auto centers = poly.get_centers();

    bool foundAny = false;  // alert if don't find any
    Thread(poly.n_faces())
        .run([&flags, &poly, &foundAny, inset_dist, popout_dist, &centers, n,
              &normals](int t, int i) {
          Flag &flag = flags[t];
          auto f = poly.faces[i];
          auto v1 = f.back();

          for (auto &v : f) {
            auto v2 = v;

            if (f.size() == n || n == 0) {
              foundAny = true;

              flag.add_vertex(i4('f', i, v),
                              tween(poly.vertexes[v], centers[i], inset_dist) +
                                  (popout_dist * normals[i]));

              flag.add_face({i4(v1), i4(v2), i4('f', i, v2), i4('f', i, v1)});
              // new inset, extruded face
              flag.add_face(i4('ex', i), i4('f', i, v1), i4('f', i, v2));
            } else {
              flag.add_face(i4(i), i4(v1), i4(v2));  // same old flag, if non-n
            }

            v1 = v2;  // current becomes previous
          }
        });

    if (!foundAny) printf("No %d - fold components were found.", n);

    flag.combine(flags);

    if (!flag.valid)
      return poly;
    else
      return Polyhedron("n" + (n ? str(n) : "") + poly.name, flag.vertexes,
                        flag.faces)
          .normalize();
  }

  // extrudeN
  // ------------------------------------------------------
  // for compatibility with older operator spec
  static Polyhedron extrudeN(Polyhedron &poly, int n = 0) {
    auto newpoly = insetN(poly, n, 0.0, 0.1);
    newpoly.name = "x" + (n ? str(n) : "") + poly.name;
    return newpoly;
  }

  // loft
  // ------------------------------------------------------
  static Polyhedron loft(Polyhedron &poly, int n = 0, float alpha = 0) {
    auto newpoly = insetN(poly, n, alpha, 0.0);
    newpoly.name = "l" + (n ? str(n) : "") + poly.name;
    return newpoly;
  }

  // Hollow (skeletonize)
  // ------------------------------------------------------
 
  static Polyhedron hollow(Polyhedron &poly, float inset_dist = 0.2,
                           float thickness = 0.1) {
    Flag flag(poly.vertexes);
    vector<Flag> flags(Thread::getnthreads());

    auto avgnormals = poly.avg_normals();
    auto centers = poly.get_centers();

    Thread(poly.n_faces())
        .run([&flags, &poly, inset_dist, thickness, &centers, &avgnormals](int t,
                                                                        int i) {
          Flag &flag = flags[t];
          auto v1 = poly.faces[i].back();

          for (auto &v2 : poly.faces[i]) {
            // new inset vertex for every vert in face
            flag.add_vertex(i4('fin', i, 'v', v2),
                            tween(poly.vertexes[v2], centers[i], inset_dist));
            flag.add_vertex(i4('fdwn', i, 'v', v2),
                            tween(poly.vertexes[v2], centers[i], inset_dist) -
                                (thickness * avgnormals[i]));

            flag.add_face(
                {i4(v1), i4(v2), i4('fin', i, 'v', v2), i4('fin', i, 'v', v1)});

            flag.add_face({i4('fin', i, 'v', v1), i4('fin', i, 'v', v2),
                           i4('fdwn', i, 'v', v2), i4('fdwn', i, 'v', v1)});
            v1 = v2;  // current becomes previous
          }
        });

    flag.combine(flags);

    if (!flag.valid)
      return poly;
    else
      return Polyhedron("H" + poly.name, flag.vertexes, flag.faces).normalize();
  }

  // Perspectiva 1
  // ------------------------------------------------------------------------------------------
  // an operation reverse-engineered from Perspectiva Corporum Regularium
  static Polyhedron perspectiva1(Polyhedron &poly) {
    auto centers = poly.get_centers();  // calculate face centers

    Flag flag;
    vector<Flag> flags(Thread::getnthreads());
    flag.set_vertexes(poly.vertexes);

    // iterate over triplets of faces v1,v2,v3
    Thread(poly.n_faces()).run([&flags, &poly, &centers](int t, int i) {
      Flag &flag = flags[t];
      auto &f = poly.faces[i];
      auto flen = f.size();
      auto v1 = f[flen - 2], v2 = f[flen - 1];
      auto vert1 = poly.vertexes[v1], vert2 = poly.vertexes[v2];

      vector<Int4> vi4;
      for (auto &v3 : f) {
        auto vert3 = poly.vertexes[v3];
        auto v12 = i4(v1, v2);  // names for "oriented" midpoints
        auto v21 = i4(v2, v1);
        auto v23 = i4(v2, v3);

        // on each Nface, N new points inset from edge midpoints towards
        // center = "stellated" points
        flag.add_vertex(v12, midpoint(midpoint(vert1, vert2), centers[i]));

        // inset Nface made of new, stellated points
        vi4.push_back(v12);

        // new tri face constituting the remainder of the stellated Nface
        flag.add_face({v23, v12, i4(v2)});

        // one of the two new triangles replacing old edge between v1->v2
        flag.add_face({i4(v1), v21, v12});

        v1 = v2;
        v2 = v3;  //  [v1, v2,[v2, v3];  // current becomes previous

        vert1 = vert2;
        vert2 = vert3;  // [vert1, vert2,[vert2, vert3];
      }
      flag.add_face(vi4);
    });

    flag.combine(flags);

    if (!flag.valid)
      return poly;
    else
      return Polyhedron("P" + poly.name, flag.vertexes, flag.faces).normalize();
  }

  //===================================================================================================
  // Goldberg-Coxeter Operators  (in progress...)
  //===================================================================================================

  // Triangular Subdivision Operator
  // ----------------------------------------------------------------------------------------------
  // limited version of the Goldberg-Coxeter u_n operator for triangular
  // meshes We subdivide manually here, instead of using the usual flag
  // machinery.
  static Polyhedron trisub(Polyhedron &poly, int n = 2) {
    for (size_t fn = 0; fn < poly.n_faces();
         fn++)  // No-Op for non-triangular meshes.
      if (poly.faces[fn].size() != 3) return poly;

    // Calculate redundant set of new vertices for subdivided mesh.
    Vertexes newVs;
    map<string, int> vmap;
    int pos = 0;

    for (size_t fn = 0; fn < poly.faces.size(); fn++) {
      auto &f = poly.faces[fn];
      auto flen = f.size();
      auto i1 = f[flen - 3], i2 = f[flen - 2],
           i3 = f[flen - 1];  // let [i1, i2, i3,f.slice(-3);
      auto v1 = poly.vertexes[i1], v2 = poly.vertexes[i2],
           v3 = poly.vertexes[i3];
      auto v21 = v2 - v1;
      auto v31 = v3 - v1;
      for (int i = 0; i <= n; i++) {
        for (int j = 0; j + i <= n; j++) {
          auto v = (v1 + v21 * (i * 1.0f / n)) + (v31 * (j * 1.0f / n));
          vmap["v" + str(fn) + "-" + str(i) + "-" + str(j)] = pos++;
          newVs.push_back(v);
        }
      }
    }

    // The above vertices are redundant along original edges,
    // we need to build an index map into a uniqueified list of them.
    // We identify vertices that are closer than a certain epsilon distance.
    float EPSILON_CLOSE = 1.0e-8f;
    Vertexes uniqVs;
    int newpos = 0;
    map<int, int> uniqmap;
    int i = 0;
    for (auto v : newVs) {
      if (uniqmap.find(i) != uniqmap.end()) continue;  // already mapped
      uniqmap[i] = newpos;
      uniqVs.push_back(v);
      for (size_t j = i + 1; j < newVs.size(); j++) {
        auto w = newVs[j];
        if (abs((v - w).len()) < EPSILON_CLOSE) uniqmap[int(j)] = newpos;
      }
      newpos++;
    }

    Faces faces;
    for (size_t fn = 0; fn < poly.n_faces(); fn++) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j + i < n; j++) {
          faces.push_back(Face{
              uniqmap[vmap["v" + str(fn) + "-" + str(i) + "-" + str(j)]],
              uniqmap[vmap["v" + str(fn) + "-" + str(i + 1) + "-" + str(j)]],
              uniqmap[vmap["v" + str(fn) + "-" + str(i) + "-" + str(j + 1)]]});
        }
      }
      for (auto i = 1; i < n; i++) {
        for (auto j = 0; j + i < n; j++) {
          faces.push_back(Face{
              uniqmap[vmap["v" + str(fn) + "-" + str(i) + "-" + str(j)]],
              uniqmap[vmap["v" + str(fn) + "-" + str(i) + "-" + str(j + 1)]],
              uniqmap[vmap["v" + str(fn) + "-" + str(i - 1) + "-" +
                           str(j + 1)]]});
        }
      }
    }

    // Create new polygon out of faces and unique vertices.
    return Polyhedron("u" + str(n) + poly.name, uniqVs, faces).normalize();
  }
};

