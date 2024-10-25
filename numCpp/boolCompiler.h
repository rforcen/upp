/*
        Logical expression compiler / executer
        with near c++ native performance < 2:1
*/
#pragma once

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <cstddef>

using namespace std;

template <typename T>
class BoolCompiler {
 public:
  enum {
    SNULL = 0,
    NUMBER = 1,

    IDENT_x = 3,
    PLUS = 5,
    MINUS = 6,
    MULT = 7,
    DIV = 8,
    OPAREN = 9,
    CPAREN = 10,
    POWER = 12,
    COMMA = 13,

    EQ = 14,
    NE = 15,
    LT = 16,
    LE = 17,
    GT = 18,
    GE = 19,

    AND = 20,
    OR = 21,
    NOT = 22,

    // function names
    FSIN = 90,
    FCOS = 91,
    FTAN = 92,
    FEXP = 93,
    FLOG = 94,
    FLOG10 = 95,
    FINT = 96,
    FSQRT = 97,
    FASIN = 98,
    FACOS = 99,
    FATAN = 100,
    FABS = 101,

    SPI = 103,

    PUSHC = 112,
    PUSHX = 113,

    NEG = 116
  };

  vector<string> fname = {"SIN",  "COS",  "TAN",  "EXP",  "LOG", "LOG10", "INT",
                          "SQRT", "ASIN", "ACOS", "ATAN", "ABS", "PI"};

  string s;  // expression to evaluate;
  int is = 0;

  int sym = 0;  // actual sym
  char ch = 0;  // actual ch
  T nval = 0;   // actual numerical value

  string id;  // actual id,

  vector<T> consts;

  bool err = false;

  string errMessage;
  int seed = 0;

  // compiler
  byte Code[4096 * 4];
  int pc = 0, CodeSize = 0;

 private:
  string toLower(const string& str) {
    string result = "";
    for (char ch : str) {
      result += tolower(ch);
    }
    return result;
  }

  void toUpper(string& s) {
    for (auto& c : s) c = toupper(c);
  }
  int indexOf(vector<string> ids, string id) {
    auto ix = std::find(ids.begin(), ids.end(), id);
    if (ix == ids.end())
      return -1;
    else
      return ix - ids.begin();
  }

 public:
  string getErrorMessage() { return errMessage; }

  bool compile(string expr) {
    pc = is = 0;
    id = "";
    s = expr;
    getch();
    err = false;

    getsym();
    Ce0();

    CodeSize = pc;
    return !err;
  }

  bool evaluate(T x) {
    // run time
    T stack[50];
    int sp = 0;

    for (int pc = 0; pc < CodeSize; pc++) {
      switch ((int)Code[pc]) {
        case PUSHC:
          stack[sp++] = consts[(int)Code[++pc]];
          break;
        case PUSHX:
          stack[sp++] = x;
          break;
        case EQ:
          sp--;
          stack[sp - 1] = stack[sp - 1] == stack[sp];
          break;
        case NE:
          sp--;
          stack[sp - 1] = stack[sp - 1] != stack[sp];
          break;
        case LT:
          sp--;
          stack[sp - 1] = stack[sp - 1] < stack[sp];
          break;
        case LE:
          sp--;
          stack[sp - 1] = stack[sp - 1] <= stack[sp];
          break;
        case GT:
          sp--;
          stack[sp - 1] = stack[sp - 1] > stack[sp];
          break;
        case GE:
          sp--;
          stack[sp - 1] = stack[sp - 1] >= stack[sp];
          break;
        case AND:
          sp--;
          stack[sp - 1] = stack[sp - 1] && stack[sp];
          break;
        case OR:
          sp--;
          stack[sp - 1] = stack[sp - 1] || stack[sp];
          break;
        case NOT:
          stack[sp - 1] = !stack[sp - 1];
          break;
        case PLUS:
          sp--;
          stack[sp - 1] += stack[sp];
          break;
        case MINUS:
          sp--;
          stack[sp - 1] -= stack[sp];
          break;
        case MULT:
          sp--;
          stack[sp - 1] *= stack[sp];
          break;
        case DIV:
          sp--;
          stack[sp - 1] /= stack[sp];
          break;
        case POWER:
          sp--;
          stack[sp - 1] = pow(stack[sp - 1], stack[sp]);
          break;
        case NEG:
          stack[sp - 1] = -stack[sp - 1];
          break;
        case FSIN:
          stack[sp - 1] = sin(stack[sp - 1]);
          break;
        case FCOS:
          stack[sp - 1] = cos(stack[sp - 1]);
          break;
        case FTAN:
          stack[sp - 1] = tan(stack[sp - 1]);
          break;
        case FASIN:
          stack[sp - 1] = asin(stack[sp - 1]);
          break;
        case FACOS:
          stack[sp - 1] = acos(stack[sp - 1]);
          break;
        case FATAN:
          stack[sp - 1] = atan(stack[sp - 1]);
          break;
        case FINT:
          stack[sp - 1] = abs(stack[sp - 1]);
          break;
        case FEXP:
          stack[sp - 1] = exp(stack[sp - 1]);
          break;
        case FLOG:
          stack[sp - 1] = log(stack[sp - 1]);
          break;
        case FLOG10:
          stack[sp - 1] = log10(stack[sp - 1]);
          break;
        case FSQRT:
          stack[sp - 1] = sqrt(stack[sp - 1]);
          break;

        default:
          err = true;
          pc = CodeSize;
          break;
      }
    }

    if (sp != 0)
      return stack[sp - 1] != 0;
    else
      return false;
  }

 private:
  bool islower(char c) { return (c >= 'a' && c <= 'z'); }
  bool isupper(char c) { return (c >= 'A' && c <= 'Z'); }
  bool isalpha(char c) { return (islower(c) || isupper(c)); }
  bool isdigit(char c) { return (c >= '0' && c <= '9'); }
  bool isalnum(char c) { return (isalpha(c) || isdigit(c)); }
  
  // get next char from *s
  char getch() {
    ch = (is < s.size()) ? s[is++] : 0;
    return ch;
  }
  void ungetch() { is--; }

  // get next symbol
  int getsym() {
    int i;

    sym = SNULL;
    id = "";

    // skip blanks
    while (ch != 0 && ch <= ' ') getch();
    
    // detect symbol
    if (isalpha(ch)) {  // ident
      id = "";
      for (i = 0; isalnum(ch) || ch == '_'; i++) {
        id += ch;
        getch();
      }

      toUpper(id);  // look up for 'x' or 't

      if (id == "X")
        sym = IDENT_x;
      else {
        // is a func ?
        int ix = 0;
        if ((ix = indexOf(fname, id)) != -1) {
          sym = ix + FSIN;  // first symbol offset
        } else {
          sym = 0;
          error("unknown symbol:" + id);
        }
      }
    } else {
      if (isdigit(ch)) {  // number (double) take care of dddd.ddde-dd
        for (i = 0; isdigit(ch) || ch == '.' || ch == 'e' || ch == 'E'; i++) {
          id += ch;
          getch();
        }
        sym = NUMBER;
        try {
          nval = (T)stod(id);
        } catch (const std::exception& e) {
          nval = 0;
          error("malformed number:" + id);
        }
      } else {
        switch (ch) {
          case '+':
            sym = PLUS;
            break;
          case '-':
            sym = MINUS;
            break;
          case '*':
            sym = MULT;
            break;
          case '/':
            sym = DIV;
            break;
          case '(':
            sym = OPAREN;
            break;
          case ')':
            sym = CPAREN;
            break;
          case '^':
            sym = POWER;
            break;
          case ',':
            sym = COMMA;
            break;

          case '&':
            sym = AND;
            break;
          case '|':
            sym = OR;
            break;

          case '=':
            if (getch() == '=')
              sym = EQ;
            else
              sym = SNULL;  // error
            break;
          case '!':
            if (getch() == '=')
              sym = NE;
            else {
              ungetch();
              sym = NOT;
            }
            break;
          case '<':
            if (getch() == '=')
              sym = LE;
            else {
              sym = LT;
              ungetch();
            }
            break;
          case '>':
            if (getch() == '=')
              sym = GE;
            else {
              sym = GT;
              ungetch();
            }
            break;

          case 0:
            sym = SNULL;
            break;

          default:
            sym = SNULL;
            error(string("character not recognized: ") + ch);
            break;
        }
        getch();
      }
    }
    return sym;
  }
  void error(string text) {
    errMessage = text;
    err = true;
  }
  void gen(int token, T f) {  // code Generation
    Code[pc++] = (byte)token;
    Code[pc++] = (byte)consts.size();
    consts.push_back(f);
  }
  void gen(int token, byte i) {
    Code[pc++] = (byte)token;
    Code[pc++] = (byte)i;
  }
  void gen(int token) { Code[pc++] = (byte)token; }
  void Ce0() {
    if (!err) {
      Ce00();
      do {
        auto _sym = sym;
        switch (sym) {
          case PLUS:
          case MINUS:
            getsym();
            Ce00();
            gen(_sym);
            break;
        }
      } while (sym == PLUS || sym == MINUS);
    }
  }
  void Ce00() {
    if (!err) {
      Ce1();
      do {
        auto _sym = sym;
        switch (sym) {
          case AND:
          case OR:
            getsym();
            Ce1();
            gen(_sym);
            break;
        }
      } while (sym == AND || sym == OR);
    }
  }

  void Ce1() {
    if (!err) {
      Ce2();
      do {
        auto _sym = sym;
        switch (sym) {
          case MULT:
          case DIV:
            getsym();
            Ce2();
            gen(_sym);
            break;
        }
      } while (sym == MULT || sym == DIV);
    }
  }
  void Ce2() {
    if (!err) {
      Ce3();
      do {
        if (sym == POWER) {
          getsym();
          Ce3();
          gen(POWER);
        }
      } while (sym == POWER);
    }
  }
  void Ce3() {
    if (!err) {
      Ce5();
      do {
        auto _sym = sym;
        switch (sym) {
          case EQ:
          case NE:
          case GT:
          case GE:
          case LT:
          case LE:
            getsym();
            Ce5();
            gen(_sym);
            break;
        }
      } while (sym == EQ || sym == NE || sym == GT || sym == GE || sym == LT ||
               sym == LE);
    }
  }

  void Ce5() {
    if (!err) {
      switch (sym) {
        case OPAREN:
          getsym();
          Ce0();
          getsym();
          break;
        case NUMBER:
          gen(PUSHC, nval);
          getsym();
          break;

        case IDENT_x:
          gen(PUSHX);
          getsym();
          break;
        case MINUS:
          getsym();
          Ce5();
          gen(NEG);
          break;
        case NOT:
          getsym();
          Ce5();
          gen(NOT);
          break;
        case PLUS:
          getsym();
          Ce5();
          break;

        case FSIN:
        case FCOS:
        case FTAN:
        case FASIN:
        case FACOS:
        case FATAN:
        case FEXP:
        case FINT:
        case FABS:
        case FLOG:
        case FLOG10:
        case FSQRT: {
          int tsym = sym;
          getsym();
          Ce5();
          gen(tsym);
        } break;

        case SPI:
          getsym();
          gen(PUSHC, M_PI);
          break;

        case SNULL:
          break;

        default:
          error("unknown symbol: " + id);
          break;  // syntax error
      }
    }
  };
};

