//
//  poly_operations.hpp
//  test_polygon
//
//  Created by asd on 05/09/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#ifndef poly_operations_hpp
#define poly_operations_hpp

#include "common.hpp"
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

    Flag flag;

    auto &mp = flag.map4;
    auto &vm = flag.vi4v;

    auto normals = poly.get_normals();
    auto centers = poly.get_centers();

    bool foundAny = false;

    // create vertex map
    auto n_vertex = flag.add_vertexes(poly.vertexes);

    // create face map
    for (size_t nface = 0; nface < poly.n_faces; nface++) {

      auto &face = poly.faces[nface];
      int v1 = face.back();

      for (auto v2 : face) {

        if (face.size() == n || n == 0) {
          foundAny = true;

          vm[i4(n_vertex + nface)] =
              flag.vtx_index(centers[nface] + (normals[nface] * apexdist));

          mp[i4(nface, v1)][i4(v1)] = i4(v2);
          mp[i4(nface, v1)][i4(v2)] = i4(n_vertex + nface);
          mp[i4(nface, v1)][i4(n_vertex + nface)] = i4(v1);

        } else {
          mp[i4(nface)][i4(v1)] = i4(v2);
        }

        v1 = v2; // current becomes previous
      }
    }

    if (!foundAny)
      printf("kisN: No %d-fold components were found.\n", n);

    return flag.to_poly("k" + (n ? str(n) : "") + poly.name);
  }

  // Ambo
  // ------------------------------------------------------
  // The best way to think of the ambo operator is as a topological "tween"
  // between a polyhedron and its dual polyhedron.  Thus the ambo of a dual
  // polyhedron is the same as the ambo of the original. Also called "Rectify".
  //

  static Polyhedron ambo(Polyhedron &poly) {
    Flag flag;

    auto &mp = flag.map4;
    auto &vm = flag.vi4v;

    for (int nface = 0; nface < poly.n_faces; nface++) {
      auto &face = poly.faces[nface];
      auto flen = face.size();

      auto v1 = face[flen - 2], v2 = face[flen - 1]; //  [v1, v2] = f.slice(-2);

      for (auto &v3 : face) {
        auto m12 = i4_min(v1, v2), m23 = i4_min(v2, v3);

        if (v1 < v2) // vertices are the midpoints of all edges of original poly
          vm[m12] =
              flag.vtx_index(midpoint(poly.vertexes[v1], poly.vertexes[v2]));

        // two new flags:
        // One whose face corresponds to the original f:
        mp[i4('orig', nface)][m12] = m23;

        // Another flag whose face  corresponds to (the truncated) v2:
        mp[i4('dual', v2)][m23] = m12;

        // shift over one
        v1 = v2;
        v2 = v3;
      }
    }

    return flag.to_poly("a" + poly.name);
  }

  // Gyro
  // -------------------------------------------------------
  // This is the dual operator to "snub", i.e dual*Gyro = Snub.  It is a bit
  // easier to implement this way.
  //
  // Snub creates at each vertex a new face, expands and twists it, and adds two
  // new triangles to replace each edge.

  static Polyhedron gyro(Polyhedron &poly) {
    Flag flag;

    auto &mp = flag.map4;
    auto &vm = flag.vi4v;

    flag.add_vertexes(poly.vertexes);

    Vertexes centers =
        poly.get_centers(); // new vertices in center of each face

    for (size_t i = 0; i < poly.n_faces; i++)
      vm[i4('cntr', i)] = flag.vtx_index(centers[i]);

    for (size_t i = 0; i < poly.n_faces; i++) {
      auto &f = poly.faces[i];
      auto flen = f.size();

      auto v1 = f[flen - 2], v2 = f[flen - 1]; //  [v1, v2] = f.slice(-2);

      for (size_t j = 0; j < flen; j++) {
        auto sv1 = str(v1), sv2 = str(v2), si = str(i);
        auto v = f[j];
        auto v3 = v;
        auto sv3 = str(v3);

        vm[i4(v1, v2)] = flag.vtx_index(
            oneThird(poly.vertexes[v1], poly.vertexes[v2])); // new v in face

        auto fname = i4(i, 'f', v1); // si + "f" + sv1;

        mp[fname][i4('cntr', i)] = i4(v1, v2); // five new flags
        mp[fname][i4(v1, v2)] = i4(v2, v1);
        mp[fname][i4(v2, v1)] = i4(v2);
        mp[fname][i4(v2)] = i4(v2, v3);
        mp[fname][i4(v2, v3)] = i4('cntr', i);

        // shift over one
        v1 = v2;
        v2 = v3;
      }
    }

    return flag.to_poly("g" + poly.name);
  }

  // Propellor
  // ------------------------------------------------------
  // builds a new 'skew face' by making new points along edges, 1/3rd the
  // distance from v1->v2, then connecting these into a new inset face.  This
  // breaks rotational symmetry about the faces, whirling them into gyres

  static Polyhedron propellor(Polyhedron &poly) {

    Flag flag;
    auto &mp = flag.map4;
    auto &vm = flag.vi4v;

    flag.add_vertexes(poly.vertexes);

    Vertexes centers =
        poly.get_centers(); // new vertices in center of each face

    for (size_t i = 0; i < poly.n_faces; i++) {

      auto &f = poly.faces[i];
      auto flen = f.size();
      auto v1 = f[flen - 2], v2 = f[flen - 1]; //  [v1, v2] = f.slice(-2);

      for (auto &v : f) {
        auto v3 = v;
        auto sv1 = str(v1), sv2 = str(v2), si = str(i), sv3 = str(v3);

        vm[i4(v1, v2)] = flag.vtx_index(
            oneThird(poly.vertexes[v1],
                     poly.vertexes[v2])); // new v in face, 1/3rd along edge

        auto fname = i4(i, 'f', v2);
        mp[i4(i)][i4(v1, v2)] = i4(v2, v3); // five new flags

        mp[fname][i4(v1, v2)] = i4(v2, v1);
        mp[fname][i4(v2, v1)] = i4(v2);
        mp[fname][i4(v2)] = i4(v2, v3);
        mp[fname][i4(v2, v3)] = i4(v1, v2);

        // shift over one
        v1 = v2;
        v2 = v3;
      }
    }

    return flag.to_poly("p" + poly.name);
  }

  // Reflection
  // ------------------------------------------------------
  // geometric reflection through origin
  static Polyhedron reflect(Polyhedron &poly) {

    // reflect each point through origin
    for (size_t i = 0; i < poly.n_vertex; i++)
      poly.vertexes[i] = -poly.vertexes[i];

    // repair clockwise-ness of faces
    for (size_t i = 0; i < poly.n_faces; i++)
      reverse(poly.faces[i].begin(), poly.faces[i].end());

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
    Flag flag;
    auto &mp = flag.map4;
    auto &vm = flag.vi4v;

    map<Flag::Int4, int> face; // make table of face as fn of edge

    for (int i = 0; i < poly.n_faces; i++) {
      auto &f = poly.faces[i];
      auto v1 = f.back(); // previous vertex index
      for (auto v2 : f) {
        face[i4(v1, v2)] = i;
        v1 = v2; // current becomes previous
      }
    }

    auto centers = poly.get_centers();
    for (int i = 0; i < poly.n_faces; i++)
      vm[i4(i)] = flag.vtx_index(centers[i]);

    for (int i = 0; i < poly.n_faces; i++) {
      auto &f = poly.faces[i];
      auto v1 = f.back(); // previous vertex
      for (auto &v2 : f) {
        mp[i4(v1)][i4(face[i4(v2, v1)])] = i4(i);
        v1 = v2; // current becomes previous
      }
    }

    auto &pn = poly.name;
    auto name = (pn[0] != 'd') ? "d" + pn : pn.substr(1, string::npos);
    return flag.to_poly(name); // build topological dual from flags
  }

  // Chamfer
  // ----------------------------------------------------
  // A truncation along a polyhedron's edges.
  // Chamfering or edge-truncation is similar to expansion, moving faces apart
  // and outward, but also maintains the original vertices. Adds a new hexagonal
  // face in place of each original edge. A polyhedron with e edges will have a
  // chamfered form containing 2e new vertices, 3e new edges, and e new
  // hexagonal faces. -- Wikipedia See also
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

  static Polyhedron chamfer(Polyhedron &poly, float dist = 0.5) {
    Flag flag;
    auto &mp = flag.map4;
    auto &vm = flag.vi4v;

    auto normals = poly.get_normals();

    // For each face f in the original poly
    for (size_t i = 0; i < poly.n_faces; i++) {
      auto &f = poly.faces[i];
      auto v1 = f.back();
      auto v1new = i4(i, v1);

      for (auto &v2 : f) {
        // TODO: figure out what distances will give us a planar hex face.
        // Move each old vertex further from the origin.
        vm[i4(v2)] = flag.vtx_index((1.0f + dist) * poly.vertexes[v2]);
        // Add a new vertex, moved parallel to normal.
        auto v2new = i4(i, v2);

        vm[v2new] =
            flag.vtx_index(poly.vertexes[v2] + (dist * 1.5f * normals[i]));

        // Four new flags:
        // One whose face corresponds to the original face:
        mp[i4('orig', i)][v1new] = v2new;

        // And three for the edges of the new hexagon:
        auto facename = v1 < v2 ? i4('hex', v1, v2) : i4('hex', v2, v1);
        mp[facename][i4(v2)] = v2new;
        mp[facename][v2new] = v1new;
        mp[facename][v1new] = i4(v1);

        v1 = v2;
        v1new = v2new;
      }
    }

    return flag.to_poly("c" + poly.name);
  }

  // Whirl
  // -------------------------------------------------------
  // Gyro followed by truncation of vertices centered on original faces.
  // This create 2 new hexagons for every original edge.
  // (https://en.wikipedia.org/wiki/Conway_polyhedron_notation#Operations_on_polyhedra)
  //
  // Possible extension: take a parameter n that means only whirl n-sided faces.
  // If we do that, the flags marked #* below will need to have their other
  // sides filled in one way or another, depending on whether the adjacent face
  // is whirled or not.

  static Polyhedron whirl(Polyhedron &poly, int n = 0) {
    (void)n;

    Flag flag;
    auto &mp = flag.map4;
    auto &vm = flag.vi4v;

    flag.add_vertexes(poly.vertexes);

    // new vertices around center of each face
    auto centers = poly.get_centers();

    for (size_t i = 0; i < poly.n_faces; i++) {
      auto &f = poly.faces[i];
      auto flen = f.size();
      auto v1 = f[flen - 2], v2 = f[flen - 1]; //  [v1, v2] = f.slice(-2);

      for (size_t j = 0; j < flen; j++) {
        auto v = f[j];
        auto v3 = v;

        // New vertex along edge
        auto v1_2 = oneThird(poly.vertexes[v1], poly.vertexes[v2]);
        vm[i4(v1, v2)] = flag.vtx_index(v1_2);
        // New vertices near center of face

        auto cv1name = i4('cntr', i, v1);
        auto cv2name = i4('cntr', i, v2);

        vm[cv1name] = flag.vtx_index(unit(oneThird(centers[i], v1_2)));

        auto fname = i4(i, 'f', v1);
        // New hexagon for each original edge
        mp[fname][cv1name] = i4(v1, v2);
        mp[fname][i4(v1, v2)] = i4(v2, v1);
        mp[fname][i4(v2, v1)] = i4(v2);
        mp[fname][i4(v2)] = i4(v2, v3);
        mp[fname][i4(v2, v3)] = cv2name;
        mp[fname][cv2name] = cv1name;

        // New face in center of each old face
        mp[i4('c', i)][cv1name] = cv2name;

        v1 = v2;
        v2 = v3;
      }
    } // shift over one

    return flag.to_poly("w" + poly.name);
  }

  // Quinto
  // -------------------------------------------------------
  // This creates a pentagon for every point in the original face, as well as
  // one new inset face.
  static Polyhedron quinto(Polyhedron &poly) {
    Flag flag;
    auto &vm = flag.vi4v;

    auto centers = poly.get_centers();

    // create vertex map
    for (int nface = 0; nface < poly.n_faces; nface++) {
      auto &f = poly.faces[nface];
      auto flen = f.size();
      auto centroid = centers[nface];

      // walk over face vertex-triplets
      auto v1 = f[flen - 2], v2 = f[flen - 1]; //  [v1, v2] = f.slice(-2);

      for (auto v3 : f) {
        // for each face-corner, we make two new points:
        Vertex midpt = midpoint(poly.vertexes[v1], poly.vertexes[v2]),
               innerpt = midpoint(midpt, centroid);

        vm[i4_min(v1, v2)] = flag.vtx_index(midpt);
        vm[i4_min(nface, v1, v2)] = flag.vtx_index(innerpt);

        // and add the old corner-vertex
        vm[i4(v2)] = flag.vtx_index(poly.vertexes[v2]);

        // shift over one
        v1 = v2;
        v2 = v3;
      }
    }
    flag.index_vertexes();

    // For each face f in the original poly
    for (int nface = 0; nface < poly.n_faces; nface++) {
      auto &f = poly.faces[nface];
      auto flen = f.size();

      // walk over face vertex-triplets
      auto v1 = f[flen - 2], v2 = f[flen - 1]; //  [v1, v2] = f.slice(-2);

      vector<Flag::Int4> vi4;

      for (auto v3 : f) {
        auto t12 = i4_min(v1, v2), ti12 = i4_min(nface, v1, v2),
             t23 = i4_min(v2, v3), ti23 = i4_min(nface, v2, v3), iv2 = i4(v2);

        // pentagon for each vertex in original face
        flag.add_face({ti12, t12, iv2, t23, ti23});
 
        // inner rotated face of same vertex-number as original
        vi4.push_back(ti12);

        // shift over one
        v1 = v2;
        v2 = v3;
      }

      flag.add_face(vi4);
    }

    return Polyhedron("q" + poly.name, flag.vertexes, flag.faces);
  }

  // inset / extrude / "Loft" operator
  // ------------------------------------------------------
  static Polyhedron insetN(Polyhedron &poly, int n = 0, float inset_dist = 0.5,
                           float popout_dist = -0.2f) {

    Flag flag;
    auto &mp = flag.map4;
    auto &vm = flag.vi4v;

    flag.add_vertexes(poly.vertexes);

    auto normals = poly.get_normals();
    auto centers = poly.get_centers();

    for (size_t i = 0; i < poly.n_faces;
         i++) { // new inset vertex for every vert in face
      auto &f = poly.faces[i];
      if (f.size() == n || n == 0) {
        for (auto &v : f) {
          vm[i4('f', i, v)] =
              flag.vtx_index(tween(poly.vertexes[v], centers[i], inset_dist) +
                             (popout_dist * normals[i]));
        }
      }
    }

    bool foundAny = false; // alert if don't find any
    for (size_t i = 0; i < poly.n_faces; i++) {
      auto f = poly.faces[i];
      auto v1 = f.back();

      for (auto &v : f) {
        auto v2 = v;

        if (f.size() == n || n == 0) {
          foundAny = true;
          auto fname = i4(i, v1);
          mp[fname][i4(v1)] = i4(v2);
          mp[fname][i4(v2)] = i4('f', i, v2);
          mp[fname][i4('f', i, v2)] = i4('f', i, v1);
          mp[fname][i4('f', i, v1)] = i4(v1);
          // new inset, extruded face
          mp[i4('ex', i)][i4('f', i, v1)] = i4('f', i, v2);
        } else {
          mp[i4(i)][i4(v1)] = i4(v2);
          // same old flag, if non-n
        }

        v1 = v2;
      }
    } // current becomes previous

    if (!foundAny)
      printf("No %d - fold components were found.", n);

    return flag.to_poly("n" + (n ? str(n) : "") + poly.name);
  }

  // extrudeN
  // ------------------------------------------------------
  // for compatibility with older operator spec
  static Polyhedron extrudeN(Polyhedron &poly, int n = 0) {
    auto newpoly = insetN(poly, n, 0.0, 0.3);
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

    Flag flag;
    auto &vm = flag.vi4v;

    auto dualnormals = PolyOperations::dual(poly).avg_normals();
    auto normals = poly.avg_normals();
    auto centers = poly.get_centers();
    auto &vertexes = poly.vertexes;

    // create vertex map v[t]=v

    for (size_t i = 0; i < poly.n_vertex;
         i++) { // each old vertex is a new vertex
      auto &p = vertexes[i];
      vm[i4('v', i)] = flag.vtx_index(p);
      vm[i4('down', i)] = flag.vtx_index(p - thickness * dualnormals[i]);
    }

    // new inset vertex for every vert in face
    for (size_t i = 0; i < poly.n_faces; i++) {
      for (auto &v : poly.faces[i]) {
        vm[i4('fin', i, 'v', v)] =
            flag.vtx_index(tween(vertexes[v], centers[i], inset_dist));
        vm[i4('fdwn', i, 'v', v)] =
            flag.vtx_index(tween(vertexes[v], centers[i], inset_dist) -
                           (thickness * normals[i]));
      }
    }
    flag.index_vertexes(); // create indexed vertex array -> flag.vertexes

    // create face map f[t][f]=t;

    for (size_t i = 0; i < poly.n_faces; i++) {
      auto &fcs = poly.faces[i];
      auto v1 = fcs.back();

      for (auto &v2 : fcs) {
        flag.add_face({i4('v', v1), i4('v', v2), i4('fin', i, 'v', v2),
                       i4('fin', i, 'v', v1)});

        flag.add_face({i4('fin', i, 'v', v1), i4('fin', i, 'v', v2),
                       i4('fdwn', i, 'v', v2), i4('fdwn', i, 'v', v1)});

        v1 = v2; // current becomes previous
      }
    }

    return Polyhedron("H" + poly.name, flag.vertexes, flag.faces);
  }

  // Perspectiva 1
  // ------------------------------------------------------------------------------------------
  // an operation reverse-engineered from Perspectiva Corporum Regularium
  static Polyhedron perspectiva1(Polyhedron &poly) {
    auto centers = poly.get_centers(); // calculate face centers

    Flag flag;
    auto &mp = flag.map4;
    auto &vm = flag.vi4v;

    flag.add_vertexes(poly.vertexes);

    // iterate over triplets of faces v1,v2,v3
    for (size_t i = 0; i < poly.n_faces; i++) {
      auto &f = poly.faces[i];
      auto flen = f.size();
      auto v1 = f[flen - 2], v2 = f[flen - 1];
      auto vert1 = poly.vertexes[v1], vert2 = poly.vertexes[v2];

      for (auto &v3 : f) {

        auto vert3 = poly.vertexes[v3];
        auto v12 = i4(v1, '~', v2); // names for "oriented" midpoints
        auto v21 = i4(v2, '~', v1);
        auto v23 = i4(v2, '~', v3);

        // on each Nface, N new points inset from edge midpoints towards center
        // = "stellated" points
        vm[v12] = flag.vtx_index(midpoint(midpoint(vert1, vert2), centers[i]));

        // inset Nface made of new, stellated points
        mp[i4('in', i)][v12] = v23;

        // new tri face constituting the remainder of the stellated Nface
        auto fname = i4('f', i, v2);
        mp[fname][v23] = v12;
        mp[fname][v12] = i4(v2);
        mp[fname][i4(v2)] = v23;

        // one of the two new triangles replacing old edge between v1->v2
        fname = i4('f', v1, '~', v2);
        mp[fname][i4(v1)] = v21;
        mp[fname][v21] = v12;
        mp[fname][v12] = i4(v1);

        v1 = v2;
        v2 = v3; //  [v1, v2] = [v2, v3];  // current becomes previous

        vert1 = vert2;
        vert2 = vert3; // [vert1, vert2] = [vert2, vert3];
      }
    }

    return flag.to_poly("P" + poly.name);
  }
  //===================================================================================================
  // Goldberg-Coxeter Operators  (in progress...)
  //===================================================================================================

  // Triangular Subdivision Operator
  // ----------------------------------------------------------------------------------------------
  // limited version of the Goldberg-Coxeter u_n operator for triangular meshes
  // We subdivide manually here, instead of using the usual flag machinery.
  static Polyhedron trisub(Polyhedron &poly, int n = 2) {

    for (size_t fn = 0; fn < poly.n_faces;
         fn++) // No-Op for non-triangular meshes.
      if (poly.faces[fn].size() != 3)
        return poly;

    // Calculate redundant set of new vertices for subdivided mesh.
    Vertexes newVs;
    map<string, int> vmap;
    int pos = 0;

    for (size_t fn = 0; fn < poly.faces.size(); fn++) {
      auto &f = poly.faces[fn];
      auto flen = f.size();
      auto i1 = f[flen - 3], i2 = f[flen - 2],
           i3 = f[flen - 1]; // let [i1, i2, i3] = f.slice(-3);
      auto v1 = poly.vertexes[i1], v2 = poly.vertexes[i2],
           v3 = poly.vertexes[i3];
      auto v21 = v2 - v1;
      auto v31 = v3 - v1;
      for (int i = 0; i <= n; i++) {
        for (int j = 0; j + i <= n; j++) {
          auto v = (v1 + (i * 1.0 / n) * v21) + (j * 1.0 / n * v31);
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
      if (uniqmap.find(i) != uniqmap.end())
        continue; // already mapped
      uniqmap[i] = newpos;
      uniqVs.push_back(v);
      for (size_t j = i + 1; j < newVs.size(); j++) {
        auto w = newVs[j];
        if (simd_distance(v, w) < EPSILON_CLOSE)
          uniqmap[int(j)] = newpos;
      }
      newpos++;
    }

    Faces faces;
    for (size_t fn = 0; fn < poly.n_faces; fn++) {
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
    return Polyhedron("u" + str(n) + poly.name, uniqVs, faces);
  }
};

#endif /* poly_operations_hpp */
