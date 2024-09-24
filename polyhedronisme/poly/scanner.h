#ifndef SCANNER_H
#define SCANNER_H

#include "common.hpp"
#include "poly_operations_mt.hpp"
#include "polyhedron.hpp"
#include "seeds.hpp"

using tokens_types = enum {
  // misc
  t_null,
  t_error,
  t_eol,
  t_comma,
  // poly
  t_Tetra,
  t_Cube,
  t_Ico,
  t_Octa,
  t_Dodeca,
  t_Prysm,
  t_Antiprism,
  t_Pyramid,
  t_Cupola,
  t_Anticupola,
  t_Johnson,
  // transformations: 'dagprPqkcwnxlH'
  t_ambo,
  t_dual,
  t_gyro,
  t_propellor,
  t_reflec,
  t_persp,
  t_quinto,
  t_kisn,
  t_chamfer,
  t_whril,
  t_insetn,
  t_extrude,
  t_loft,
  t_hollow,
  // number
  t_number
};

class Token {
 public:
  Token() : n(0), f1(0), f2(0) {}
  Token(tokens_types token) : token(token), n(0), f1(0), f2(0) {}
  Token(int n, float f1, float f2) : n(n), f1(f1), f2(f2) {}
  Token(tokens_types token, int n, float f1, float f2)
      : token(token), n(n), f1(f1), f2(f2) {}

  tokens_types token;
  int n;
  float f1, f2;
};

class Scanner {
  map<char, tokens_types> token_map = {
      // seeds
      {'T', t_Tetra},
      {'C', t_Cube},
      {'I', t_Ico},
      {'O', t_Octa},
      {'D', t_Dodeca},
      {'J', t_Johnson},
      {'P', t_Prysm},
      {'A', t_Antiprism},
      {'Y', t_Pyramid},
      {'U', t_Cupola},
      {'V', t_Anticupola},

      // trasformations
      {'a', t_ambo},
      {'d', t_dual},
      {'g', t_gyro},
      {'p', t_propellor},
      {'r', t_reflec},
      {'P', t_persp},
      {'q', t_quinto},
      {'k', t_kisn},
      {'c', t_chamfer},
      {'w', t_whril},
      {'n', t_insetn},
      {'x', t_extrude},
      {'l', t_loft},
      {'H', t_hollow}};

  map<tokens_types, Token> default_values = {
      {t_kisn, {0, 0.1, 0}},    {t_chamfer, {0, 0.05, 0}},
      {t_whril, {0, 0, 0}},     {t_insetn, {0, 0.3, -0.1}},
      {t_extrude, {0, 0, 0}},   {t_loft, {0, 0, 0}},
      {t_hollow, {0, 0.2, 0.1}}};

  map<char, string> transf_equivalents = {
      {'t', "dkd"}, {'b', "ta"},  {'e', "aa"}, {'j', "dad"},
      {'o', "jj"},  {'m', "k3j"}, {'s', "dgd"}};

 public:
  bool is_poly(tokens_types k) { return k >= t_Tetra && k <= t_Johnson; }
  bool is_transfor(tokens_types k) { return k >= t_ambo && k <= t_hollow; }
  bool is_float() { return n != int(f); }
  void set_default_value(tokens_types token, Token &tok) {
    auto dv = default_values.find(token);
    if (dv != default_values.end()) {
      auto t = dv->second;
      tok.n = t.n;
      tok.f1 = t.f1;
      tok.f2 = t.f2;
    }
  }

  string replace_equivalence(string s) {
    string sr = s;
    bool replaced;
    do {
      replaced = false;
      auto sw = sr;
      sr = "";
      for (auto &c : sw) {
        auto fc = transf_equivalents.find(c);
        if (fc != transf_equivalents.end()) {
          sr += fc->second;
          replaced = true;
        } else
          sr += c;
      }
    } while (replaced);
    return sr;
  }

  tokens_types get_token() {
    if (i >= l) return t_eol;

    char ch = s[i];

    if (token_map.find(ch) != token_map.end()) {
      i++;
      return token_map[ch];
    }

    if (isnumber(ch) || ch == '-' || ch == '.') {
      string num;

      do {
        num += ch;
        ch = s[++i];
      } while (
          (isnumber(ch) || ch == '.' || ch == 'e' || ch == 'E' || ch == '-') &&
          i <= l);
      n = atoi(num.c_str());
      f = atof(num.c_str());
      return t_number;
    }

    if (ch == ',') {
      i++;
      return t_comma;
    }

    i++;
    return t_error;
  }

  void init(string s) {
    this->s = s;
    i = 0;
    l = s.size();
  }

  Polyhedron scan(string s) {
    vector<Token> vk;

    init(replace_equivalence(s));

    auto token = get_token();
    do {
      if (is_transfor(token)) {
        Token t(token);

        set_default_value(token, t);

        if ((token = get_token()) == t_number) {
          bool ff = is_float();  // first float

          if (ff) {
            t.n = 0;
            t.f1 = f;
          } else {
            t.n = n;
          }

          if (((token = get_token()) == t_comma) && (get_token() == t_number)) {
            if (ff) {
              t.f2 = f;
              token = get_token();
            } else {
              t.f1 = f;
              if (((token = get_token()) == t_comma) &&
                  (get_token() == t_number)) {
                t.f2 = f;
                token = get_token();
              }
            }
          }
        }

        vk.push_back(t);

      } else if (is_poly(token)) {
        Token tk(token, 4, 0, 0);
        if ((token = get_token()) == t_number) {
          tk.n = n;
          if (((token = get_token()) == t_comma) && (get_token() == t_number)) {
            tk.f1 = f;
            if (((token = get_token()) == t_comma) &&
                (get_token() == t_number)) {
              tk.f2 = f;
              token = get_token();
            }
          }
        }
        vk.push_back(tk);
      }
    } while (token != t_eol);

    reverse(vk.begin(), vk.end());
    return vk_to_poly(vk);
  }

  Polyhedron vk_to_poly(vector<Token> &vk) {
    Polyhedron p;

    for (auto &v : vk) {
      switch (v.token) {
        // seeds
        case t_Tetra:
          p = Seeds::tetrahedron();
          break;
        case t_Cube:
          p = Seeds::cube();
          break;
        case t_Ico:
          p = Seeds::icosahedron();
          break;
        case t_Octa:
          p = Seeds::octahedron();
          break;
        case t_Dodeca:
          p = Seeds::dodecahedron();
          break;

        case t_Prysm:
          p = Seeds::prism(v.n);
          break;
        case t_Antiprism:
          p = Seeds::antiprism(v.n);
          break;
        case t_Cupola:
          p = Seeds::cupola(v.n, v.f1, v.f2);
          break;
        case t_Anticupola:
          p = Seeds::anticupola(v.n, v.f1, v.f2);
          break;
        case t_Johnson:
          p = Seeds::johnson(v.n);
          break;

          // transformations
        case t_kisn:
          p = PolyOperations::kisN(p, v.n, v.f1);
          break;
        case t_ambo:
          p = PolyOperations::ambo(p);
          break;
        case t_gyro:
          p = PolyOperations::kisN(p);
          break;
        case t_reflec:
          p = PolyOperations::reflect(p);
          break;
        case t_dual:
          p = PolyOperations::dual(p);
          break;
        case t_chamfer:
          p = PolyOperations::chamfer(p, v.f1);
          break;
        case t_whril:
          p = PolyOperations::whirl(p, v.n);
          break;
        case t_quinto:
          p = PolyOperations::quinto(p);
          break;
        case t_insetn:
          p = PolyOperations::insetN(p, v.n, v.f1, v.f2);
          break;
        case t_extrude:
          p = PolyOperations::extrudeN(p, v.n);
          break;
        case t_loft:
          p = PolyOperations::loft(p, v.n, v.f1);
          break;
        case t_hollow:
          p = PolyOperations::hollow(p, v.f1, v.f2);
          break;
        case t_persp:
          p = PolyOperations::perspectiva1(p);
          break;
        default:
          break;
      }
    }

    p.recalc();
    return p;
  }

 private:
  string s;
  int i = 0, l = 0;
  int n;
  float f;
};

#endif  // SCANNER_H
