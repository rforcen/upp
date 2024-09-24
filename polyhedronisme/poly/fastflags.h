#pragma once

#include "common.hpp"
#include "polyhedron.hpp"

class Flag {
 public:
  Vertexes vertexes;
  Faces faces;

  vector<I4Vix> v;
  vector<MapIndex> m;  // m[i4][i4]=i4 -> m[]<<i4,i4,i4
  vector<vector<Int4>> fcs;

  int v_index = 0;  // index of last added vertex (add_vertex)
  bool valid = true;

  Flag() = default;
  Flag(vector<Flag> &flags) {  // consolidate flags -> flag
    combine(flags);
  }

  Flag(Vertexes &vertexes) {  // consolidate flags -> flag
    set_vertexes(vertexes);
  }

  void combine(vector<Flag> &flags) {  // combine flags(v,m,fcs) -> flag

    // totalize v, fcs -> calc offsets
    int vs = v.size(), v_tot = vs, f_tot = 0,
        m_tot = 0;  // v may contain vertex
    vector<int> v_offsets{0}, fcs_offsets{0}, m_offsets{0};

    v_offsets[0] = vs;

    for (auto &f : flags) {
      v_offsets.push_back(v_tot += f.v.size());
      fcs_offsets.push_back(f_tot += f.fcs.size());
      m_offsets.push_back(m_tot += f.m.size());
    }

    // resize v,m to total required space
    v.resize(v_tot + v.size());
    m.resize(m_tot);

    // copy threaded flag.v,m,fcs << flags[].v,m,fcs
    Thread(flags.size())
        .run([this, &flags, &v_offsets,
              &m_offsets](int nflag) {  // combine (v,m) flags[] -> flag
          copy(flags[nflag].v.begin(), flags[nflag].v.end(),
               &v[v_offsets[nflag]]);
          copy(flags[nflag].m.begin(), flags[nflag].m.end(),
               &m[m_offsets[nflag]]);
        });

    index_vertexes();  // numerate 'v', v->vertexes

    if (process_m()) {  // faces << m, is valid?
      // faces << fcs
      int fs_m = faces.size();
      faces.resize(f_tot + fs_m);

      Thread(flags.size()).run([this, &flags, &fcs_offsets, fs_m](int t) {
        int offset = fcs_offsets[t] + fs_m;
        for (auto &fc : flags[t].fcs) {
          Face face;
          for (auto &vix : fc) face.push_back(find_vertex_index(vix));
          faces[offset++] = face;
        }
      });
    }
  }

  inline int add_vertex(Vertex v) {
    vertexes.push_back(v);
    return vertexes.size() - 1;
  }

  void sort_unique_v() {  // fastest solution
    sort(v.begin(), v.end(), [](I4Vix &a, I4Vix &b) -> bool { return a < b; });

    const auto &it = unique(v.begin(), v.end(),
                            [](I4Vix &a, I4Vix &b) -> bool { return a == b; });
    v.resize(std::distance(v.begin(), it));
  }

  void index_vertexes() {  // v, numerate vertexes index & create vertexes[]
    sort_unique_v();
    vertexes.resize(v.size());  // numerate & create vertexes[]

    Thread(v.size()).run([this](int i) {
      VertexIndex &_v = v[i].vix;  // <index, vertex>
      _v.index = i;
      vertexes[i] = _v.vertex;
    });
  }

  int set_vertexes(Vertexes &vertexes) {  // v = vertexes
    v.resize(vertexes.size());
    Thread(vertexes.size()).run([this, &vertexes](int i) {
      v[i] = {to_int4(i), {0, vertexes[i]}};
    });
    return vertexes.size();
  }

  inline void add_face(Int4 i0, Int4 i1, Int4 i2) { m.push_back({i0, i1, i2}); }
  inline void add_face(vector<Int4> v) { fcs.push_back(v); }

  inline void add_vertex(Int4 ix, Vertex vtx) {  // to v
    v.push_back({ix, {v_index++, vtx}});
  }

  inline int set_vertex(int i, Int4 ix, Vertex vtx) {
    v[i] = {ix, {0, vtx}};
    return i;
  }

  inline I4Vix find_vertex(Int4 _v) {
    return *lower_bound(v.begin(), v.end(), _v, I4Vix::less);
  }

  inline int find_vertex_index(Int4 _v) {
    return lower_bound(v.begin(), v.end(), _v, I4Vix::less)->vix.index;
  }

  inline Int4 find_m(Int4 _m0, Int4 _m1) {
    return (*lower_bound(m.begin(), m.end(), MapIndex(_m0, _m1, 0))).i2;
  }

  // gen. vector of from index of face change in m
  vector<int> from_to_m() {
    vector<int> v_ft;
    Int4 c0 = m[0].i0;
    int from = 0;

    for (int i = 0; i < m.size(); i++) {
      if (m[i].i0 != c0) {
        v_ft.push_back(from);
        from = i;
        c0 = m[i].i0;
      }
    }
    v_ft.push_back(from);

    return v_ft;
  }

  void to_poly() {  // v,m -> vertexes, faces

    index_vertexes();

    faces.clear();

    process_m();
  }

  bool process_m() {  // m->faces, -> valid
    const int MaxIters = 100;

    if (!m.empty()) {
      // sort m
      sort(m.begin(), m.end(), MapIndex::less);

      auto ft = from_to_m();  // from index in m (segments of i0)
      faces = Faces(ft.size());

      Thread(ft.size()).run([this, &ft](int fti) {
        int i = ft[fti];

        auto &m0 = m[i];
        Int4 _v0 = m0.i2, _v = _v0, _m0 = m0.i0;

        // traverse _m0
        Face face;
        int cnt = 0;
        do {
          face.push_back(find_vertex_index(_v));
          _v = find_m(_m0, _v);
        } while (_v != _v0 || cnt++ > MaxIters);

        faces[fti] = face;

        valid = (cnt >= MaxIters);
      });
    }
    
    return valid;
  }

  struct Int4int {  // face map
    Int4 _i4;
    int i;
    bool operator<(const Int4int &o) const { return _i4 < o._i4; }
    bool operator<(const Int4 &o) const { return _i4 < o; }
    bool operator()(const Int4int &a, const Int4 &b) const { return a._i4 < b; }
    static Int4 find(const vector<Int4int> &iv, Int4 k) {
      return i4(lower_bound(iv.begin(), iv.end(), k)->i);
    }
    static void sort(vector<Int4int> &iv) { std::sort(iv.begin(), iv.end()); }
  };

  static vector<Int4int> gen_face_map(Polyhedron &poly) {
    // make table of face as fn of edge

    vector<Int4int> face_map;

    for (int i = 0; i < poly.n_faces(); i++) {
      auto &f = poly.faces[i];
      auto v1 = f.back();  // previous vertex index
      for (auto v2 : f) {
        face_map.push_back({i4(v1, v2), i});
        v1 = v2;  // current becomes previous
      }
    }
    Int4int::sort(face_map);

    return face_map;
  }
};

