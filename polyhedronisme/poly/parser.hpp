//
//  parser.hpp
//  test_polygon
//
//  Created by asd on 09/09/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#ifndef parser_hpp
#define parser_hpp

#include "common.hpp"
#include "poly_operations_mt.hpp"
#include "polyhedron.hpp"
#include "seeds.hpp"
#include <ctype.h>

class Parser {
 public:
  Parser() {}

  static Polyhedron randomPoly() { return parse(randomTransform()); }

  static string randomTransform() {  // "dagprPqkcwnxlH"
    char const *_tr = "dagprPqkcwnxlH", *tr = "darqkc", *pl = "TCIODPAYUVJ";
    string str;

    for (int i = 0; i < 5 + (rand() % 10); i++)
      str += tr[rand() % strlen(tr)];  // transformations

    string p(1, pl[rand() % strlen(pl)]);  // poly
    str += (p[0] == 'J') ? p + to_string(rand() % johnsons.size() + 1) : p;

    return str;
  }

  static Polyhedron parse(string s, int vertexLimit=10000) {  // ttttBN
    Polyhedron p;
    int n = 0;
    size_t slen = s.length(), i = 0;
    string sd;

    reverse(s.begin(), s.end());  // NBtttt

    for (i = 0; isdigit(s[i]); i++) sd += s[i];  // N

    if (!sd.empty()) {
      reverse(sd.begin(), sd.end());
      try {
        n = std::stoi(sd);
      } catch (std::invalid_argument) {
        n = -1;
      }
    }

    switch (s[i]) {  //  base poly
      case 'T':
        p = Seeds::tetrahedron();
        break;
      case 'C':
        p = Seeds::cube();
        break;
      case 'O':
        p = Seeds::octahedron();
        break;
      case 'I':
        p = Seeds::icosahedron();
        break;
      case 'D':
        p = Seeds::dodecahedron();
        break;
      case 'P':
        p = n == 0 ? Seeds::prism() : Seeds::prism(n);
        break;
      case 'A':
        p = n == 0 ? Seeds::antiprism() : Seeds::antiprism(n);
        break;
      case 'Y':
        p = n == 0 ? Seeds::pyramid() : Seeds::pyramid(n);
        break;
      case 'U':
        p = n == 0 ? Seeds::cupola() : Seeds::cupola(n);
        break;
      case 'V':
        p = n == 0 ? Seeds::anticupola() : Seeds::anticupola(n);
        break;
      case 'J':
        p = Seeds::johnson(n % johnsons.size());
        break;
      default:
        return p;  // wrong base
    }

    for (i++; i < slen; i++) {  // transformations: dagprPqkcwnxlH
      switch (s[i]) {
        case 'd':
          p = PolyOperations::dual(p);
          break;
        case 'a':
          p = PolyOperations::ambo(p);
          break;
        case 'g':
          p = PolyOperations::gyro(p);
          break;
        case 'p':
          p = PolyOperations::propellor(p);
          break;
        case 'r':
          p = PolyOperations::reflect(p);
          break;
        case 'P':
          p = PolyOperations::perspectiva1(p);
          break;
        case 'q':
          p = PolyOperations::quinto(p);
          break;

        case 'k':
          p = PolyOperations::kisN(p);
          break;  // parameters
        case 'c':
          p = PolyOperations::chamfer(p);
          break;
        case 'w':
          p = PolyOperations::whirl(p);
          break;
        case 'n':
          p = PolyOperations::insetN(p);
          break;
        case 'x':
          p = PolyOperations::extrudeN(p);
          break;
        case 'l':
          p = PolyOperations::loft(p);
          break;
        case 'H':
          p = PolyOperations::hollow(p);
          break;

        default:
          break;
      }
      if (p.n_vertexes() > vertexLimit) break; 
    }

    p.scale_vertexes();
    return p.recalc();
  }
};

#endif /* parser_hpp */
