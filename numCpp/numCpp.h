#pragma once

#include <stdlib.h>
#include <vector>
#include <numeric>
#include <functional>
#include <assert.h>
#include <memory.h>
#include <algorithm>
#include <cstdarg>
#include <random>
#include <barrier>
#include <regex>
#include <fstream>

#include "Thread.h"
#include "compiler.h"

using namespace std;

typedef vector<vector<int>> VVInt;
typedef vector<int> VInt;

template <typename T>  // vector concat '+' operator
vector<T> operator+(vector<T> a, vector<T> b) {
  vector<T> result;
  result.reserve(a.size() + b.size());
  result.insert(result.end(), a.begin(), a.end());
  result.insert(result.end(), b.begin(), b.end());
  return result;
}

/////////////////////

template <typename T>
class Vec : public vector<T> {
 public:
  T prod() {
    T p = 1;
    for (auto &v : *this) p *= v;
  }
  void sort() {
    qsort(this->begin(), this->end(), [=](const void *a, const void *b) {
      return *(T)a < *(T *)b ? -1 : 1;
    });
  }
  bool operator==(Vec &o) {
    return equal(this->begin(), this->end(), o.begin());
  }
  T ix(int i) { return this->data[this->size() - i - 1]; }
};

template <typename T = double>
class NC {
 public:
  vector<T> data;
  VInt dims, mlt;
  int ndims = 0;
  int dim0 = 0;
  int size = 0, sizeBytes = 0;

  typedef struct _HeaderNPY {
    _HeaderNPY(uint16_t headLen) : headLen(headLen) {}
    _HeaderNPY() {}

    char id[6] = {'\x93', 'N', 'U', 'M', 'P', 'Y'};
    char minVer = 1, maxVer = 0;
    uint16_t headLen;
  } HeaderNPY;

 public:
  NC() {}

  NC(const NC &o) {
    data.clear();
    dims.clear();
    mlt.clear();

    data.assign(o.data.begin(), o.data.end());
    dims.assign(o.dims.begin(), o.dims.end());
    mlt.assign(o.mlt.begin(), o.mlt.end());

    ndims = o.ndims;
    dim0 = (ndims > 0) ? dim(0) : 0;
    size = o.size;
    sizeBytes = o.sizeBytes;
  }

  template <typename U>  // construct from another type
  NC(NC<U> &o) {
    data.clear();
    dims.clear();
    mlt.clear();

    dims.assign(o.dims.begin(), o.dims.end());
    mlt.assign(o.mlt.begin(), o.mlt.end());

    ndims = o.ndims;
    dim0 = (ndims > 0) ? dim(0) : 0;
    size = o.size;
    sizeBytes = o.sizeBytes;

    data.resize(size);  // alloc & copy data
    o.copyTo(*this);
  }

  template <typename... Args>
  NC(Args... args) {
    dims = VInt{args...};
    __recalc();
  }

 private:
  void __recalc() {
    ndims = dims.size();
    dim0 = (ndims > 0) ? dim(0) : 0;

    size = accumulate(dims.begin(), dims.end(), 1, multiplies<int>());
    sizeBytes = size * sizeof(T);
    data.resize(size);
    mlt.resize(ndims);
    for (int i = 0, nn = 1; i < ndims; i++) {
      mlt[ndims - 1 - i] = nn;
      nn *= dims[ndims - 1 - i];
    }
  }

 public:
  static NC arange(int r) {
    NC nc(r);
    for (int i = 0; i < r; i++) nc[i] = i;

    return nc;
  }

  static NC identity(int n) {
    NC nc(n, n);
    for (int i = 0; i < n; i++) nc(i, i) = 1;

    return nc;
  }

  ~NC() {
    data.clear();
    dims.clear();
    mlt.clear();
  }

  template <typename... Args>
  static NC random(Args... args) {
    return NC<T>(args...).randMT();
  }

  NC rand() {
    for (auto &d : data) d = static_cast<T>(::random()) / RAND_MAX;
    return *this;
  }
  NC randMT() {
    const int _max_range = 32000;
    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<int> distribution(1, _max_range);

    MyThread(size).run([&](int t, int from, int to) {
      for (int i = from; i < to; i++)
        data[i] = static_cast<T>(distribution(generator)) / _max_range;
    });

    return *this;
  }
  NC fill(T f) {
    for (auto &d : data) d = f;
    return *this;
  }
  NC zeros() {
    for (auto &d : data) d = 0;
    return *this;
  }
  T sum() {
    T s = static_cast<T>(0);
    for (auto &d : data) s += d;
    return s;
  }
  T prod() {
    T p = static_cast<T>(1 + 1e-10);
    for (auto &d : data) p *= d;
    return p;
  }

  NC lin(T from = T(0)) {
    for (auto &d : data) d = from++;
    return *this;
  }

  void disp() {
    for (auto &d : data) cout << d << ' ';
    cout << endl;
  }

  void print(string s) {
    cout << s << '(' << dimsStr() << ") size:(" << size << "):";
    if (size == 0) return;

    print(0, 0);
    cout << endl;
  }
  void print(int level = 0, int offset = 0) {
    if (level == ndims - 1) {
      cout << "[";
      for (int i = 0; i < dims[level]; ++i) {
        cout << data[offset + i];
        if (i < dims[level] - 1) cout << ", ";
      }
      cout << "]";
    } else {
      cout << "[";
      for (int i = 0; i < dims[level]; ++i) {
        print(level + 1, offset + i * mlt[level]);
        if (i < dims[level] - 1) cout << ", ";
      }
      cout << "]";
    }
  }

  NC apply(function<T(T)> const &lambda) { return NC(*this).inApplyMT(lambda); }
  NC applyMT(function<T(T)> const &lambda) {
    return NC(*this).inApplyMT(lambda);
  }
  NC inApply(function<T(T)> const &lambda) {  // in place
    for (auto &d : data) d = lambda(d);
    return *this;
  }
  NC inApplyMT(function<T(T)> const &lambda) {  // in place
    MyThread(size).run([&](int t, int from, int to) {
      auto d = data.begin() + from;
      for (int i = from; i < to; i++, d++) *d = lambda(*d);
    });
    return *this;
  }

  // sort section
  NC flatSort() { return NC(*this).inFlatSort(); }

  NC sort() {
    NC res(*this);
    if (ndims <= 2)
      res = res.flatSort();
    else {
      VInt tdim(dims.begin(), dims.end() - 2);  // dims[0..-2]
      for (auto &t : combinations(tdim)) res.assign(slice(t).flatSort(), t);
    }

    return res;
  }

  NC inFlatSort() {  // ascensing
    qsort((void *)data.data(), data.size(), sizeof(T),
          [](const void *a, const void *b) {
            return (*(T *)a < *(T *)b) ? -1 : 1;
          });
    return *this;
  }
  NC inFlatSortDesc() {  // ascensing
    qsort((void *)data.data(), data.size(), sizeof(T),
          [](const void *a, const void *b) {
            return (*(T *)a < *(T *)b) ? -1 : 1;
          });
    return *this;
  }

  NC inSort() {  // in place sort
    if (ndims <= 2)
      inFlatSort();
    else {
      VInt tdim(dims.begin() + 0, dims.end() - 2);  // dims[0..-2]
      for (auto &t : combinations(tdim)) assign(slice(t).flatSort(), t);
    }
    return *this;
  }

  // data acess
  auto &getData() { return data; }
  const auto begin() { return data.begin(); }
  const auto end() { return data.end(); }

  template <typename U>  // copy from a different type U
  void copyTo(NC<U> &o) {
    transform(data.begin(), data.end(), o.begin(),
              [=](T v) { return static_cast<U>(v); });
  }

 private:
  bool areOper(NC &o) { return (size == o.size); }
  bool diffDims(const NC &o) { return (size != o.size); }
  bool sameDims(NC &o) {
    return areOper(o) && equal(dims.begin(), dims.end(), o.dims.begin());
  }

  // index
  bool eqDims(VInt &dims) { return this->ndims == dims.size(); }

  template <typename... Args>
  int calcIndex(Args... args) {
    VInt index = {args...};  // vector<int>{args...}

    if (!eqDims(index)) throw out_of_range("wrong index size");

    switch (ndims) {
      case 1:
        return index[0];
      case 2:
        return index[0] * mlt[0] + index[1];
      case 3:
        return index[0] * mlt[0] + index[1] * mlt[1] + index[2];
      case 4:
        return index[0] * mlt[0] + index[1] * mlt[1] + index[2] * mlt[2] +
               index[3];
      default: {
        int res = 0;
        for (int i = 0; i < mlt.size(); i++) res += mlt[i] * index[i];
        return res;
      }
    }
  }

  int calcIndex0(VInt index) {
    if (!eqDims(index)) throw out_of_range("wrong index size");

    switch (ndims) {
      case 1:
        return index[0];
      case 2:
        return index[0] * mlt[0] + index[1];
      case 3:
        return index[0] * mlt[0] + index[1] * mlt[1] + index[2];
      case 4:
        return index[0] * mlt[0] + index[1] * mlt[1] + index[2] * mlt[2] +
               index[3];
      default: {
        int res = 0;
        for (int i = 0; i < mlt.size(); i++) res += mlt[i] * index[i];
        return res;
      }
    }
  }

 public:
  template <typename... Args>
  T &at(Args... args) {
    int ix = calcIndex(args...);
    if (ix >= size || ix < 0) throw out_of_range("index out of range");

    return data[ix];
  }
  T &at(vector<int> args) {
    int ix = calcIndex(args);
    if (ix >= size || ix < 0) throw out_of_range("index out of range");

    return data[ix];
  }
  inline T &at(int r, int c) { return data[r * dim0 + c]; }
  inline T &at(int i) { return data[i]; }

  T &operator[](initializer_list<int> args) {
    int ix = calcIndex0(args);
    if (ix >= size || ix < 0) throw out_of_range("index out of range");

    return data[ix];
  }
  const T &operator[](initializer_list<int> args) const {
    int ix = calcIndex0(args);
    if (ix >= size || ix < 0) throw out_of_range("index out of range");

    return data[ix];
  }

  // r,c index
  inline T &operator()(int r, int c) { return data[r * dim0 + c]; }  // nc(i,j)
  inline T &operator[](int i) { return data[i]; }                    // nc(i)

  // operators (+)
  NC operator+(NC o) {
    if (diffDims(o)) throw out_of_range("incompatible dimensions");

    NC r(*this);
    for (int i = 0; i < r.size; i++) r.at(i) += o.at(i);
    return r;
  }
  void operator+=(NC &o) {
    if (diffDims(o)) throw out_of_range("incompatible dimensions");
    for (int i = 0; i < size; i++) at(i) += o(i);
  }
  NC operator+(T o) {
    NC r(*this);
    for (auto &d : r.data) d += o;
    return r;
  }
  void operator+=(T o) {
    for (auto &d : data) d += o;
  }
  // (-)
  NC operator-(NC &o) {
    assert(areOper(o) && "different dimensions");
    NC r(*this);
    for (int i = 0; i < r.size; i++) r.at(i) -= o.at(i);
    return r;
  }
  void operator-=(NC &o) {
    if (diffDims(o)) throw out_of_range("incompatible dimensions");
    for (int i = 0; i < size; i++) at(i) -= o(i);
  }
  NC operator-(T o) {
    NC r(*this);
    for (int i = 0; i < r.size; i++) r(i) -= o;
    return r;
  }
  void operator-=(T o) {
    for (auto &d : data) d -= o;
  }
  // (*)
  NC operator*(NC &o) {
    assert(areOper(o) && "different dimensions");
    NC r(*this);
    for (int i = 0; i < r.size; i++) r.at(i) *= o.at(i);
    return r;
  }
  void operator*=(NC &o) {
    if (diffDims(o)) throw out_of_range("incompatible dimensions");
    for (int i = 0; i < size; i++) at(i) *= o.at(i);
  }
  NC operator*(T o) {
    NC r(*this);
    for (int i = 0; i < r.size; i++) r(i) *= o;
    return r;
  }
  void operator*=(T o) {
    for (auto &d : data) d *= o;
  }
  // (/)
  NC operator/(NC &o) {
    if (diffDims(o)) throw out_of_range("incompatible dimensions");
    NC r(*this);
    for (int i = 0; i < r.size; i++) r(i) /= o(i);
    return r;
  }
  void operator/=(NC &o) {
    if (diffDims(o)) throw out_of_range("incompatible dimensions");
    for (int i = 0; i < size; i++) at(i) /= o(i);
  }
  NC operator/(T o) {
    NC r(*this);
    for (int i = 0; i < r.size; i++) r(i) /= o;
    return r;
  }
  void operator/=(T o) {
    for (auto &d : data) d /= o;
  }

  // logical
  bool operator==(NC &o) {
    assert(areOper(o) && "different dimensions");
    return equal(data.begin(), data.end(), o.data.begin());
  }
  bool operator!=(NC &o) { return !(o == *this); }
  bool operator>(NC &o) {
    assert(areOper(o) && "different dimensions");
    return data > o.data;
  }
  bool operator>=(NC &o) {
    assert(areOper(o) && "different dimensions");
    return data >= o.data;
  }
  bool operator<(NC &o) {
    assert(areOper(o) && "different dimensions");
    return data < o.data;
  }
  bool operator<=(NC &o) {
    assert(areOper(o) && "different dimensions");
    return data <= o.data;
    return true;
  }

  NC operator<<(NC o) {  // append 'o' to this
    return push_back(o);
  }

  NC concat(NC o) const {  // join creating a 1d concat array -> after reshape
    NC res(*this);

    res.data.insert(res.data.end(), o.data.begin(),
                    o.data.end());  // insert 'o'

    return res.reshape({res.size});  // 1d {size}
  }

  NC push_back(NC o) {                                  // in place
    if (ndims == 1 && o.ndims == 1 && o.dim0 < dim0) {  // append this << o
      data.insert(data.end(), o.data.begin(), o.data.end());
      dims[0] += o.dims[0];
      __recalc();
    } else {
      if (sameDims(o)) {  // nd + nd
        data.insert(data.end(), o.data.begin(), o.data.end());
        dims[0] *= 2;
        __recalc();
      } else {
        if (ndims - 1 == o.ndims && equal(dims.begin() + 1, dims.end(),
                                          o.dims.begin())) {  // nd + md (n=m+1)
          data.insert(data.end(), o.data.begin(), o.data.end());
          dims[0]++;
          __recalc();
        } else
          throw runtime_error("push_back: incompatible dims");
      }
    }

    return *this;
  }

  NC stackVert(NC o) {  // in place
    if (sameDims(o)) {
      data.insert(data.begin(), o.data.begin(), o.data.end());  // ins. top
      dims.insert(dims.begin(), 2);  // new dim = {2} + dims
      __recalc();
    } else if (ndims - 1 == o.ndims &&
               equal(dims.begin() + 1, dims.end(), o.dims.begin())) {
      data.insert(data.begin(), o.data.begin(), o.data.end());
      dim(0)++;
      __recalc();
    } else
      throw runtime_error("stackVert: incompatible dims");
    return *this;
  }

  NC stackHorz(NC o) {  // in place
    if (sameDims(o)) {
      data.insert(data.end(), o.data.begin(), o.data.end());  // ins. end
      dims.insert(dims.begin(), 2);  // new dim = {2} + dims
      __recalc();
    } else if (ndims - 1 == o.ndims &&
               equal(dims.begin() + 1, dims.end(), o.dims.begin())) {
      data.insert(data.end(), o.data.begin(), o.data.end());
      dim(0)++;
      __recalc();
    } else
      throw runtime_error("stackVert: incompatible dims");
    return *this;
  }

  NC filter(function<bool(T)> test) const {
    NC res(dims);
    for (int i = 0; i < size; i++)
      if (test(data[i])) res.data[i] = data[i];
    return res;
  }

  NC filter(string expr) const {  // filter by evaluated string expr.
    NC res(dims);

    Evaluator ev(expr);  // compile
    if (ev.ok) {
      for (int i = 0; i < size; i++)
        if (ev.evaluate(data[i])) res.data[i] = data[i];
    } else
      throw runtime_error(
          strcat((char *)"syntax error in expression:", (char *)expr.c_str()));
    return res;
  }

  NC match(function<bool(T)> test)
      const {  // return a 1d array of items that match test
    NC res;

    for (int i = 0; i < size; i++)
      if (test(data[i])) res.data.push_back(data[i]);

    res.dims = {(int)res.data.size()};  // 1d size
    res.__recalc();
    return res;
  }

  NC expand(VInt xdims) {  // expand to higher dimension xdims + dims

    NC o(*this), res(xdims + dims);

    for (auto &c : combinations(xdims)) res.assign(o, c);

    return res;
  }

  NC contract(int nc, function<NC(NC, VInt)>
                          fc) {  // contract by applying fc to every submatrix
                                 // fc contracts from NC.dims to VInt dims

    VInt ndim(dims), sdim;  // remove nc items from nd.begin
    ndim.erase(ndim.begin(), ndim.begin() + nc);
    sdim.insert(sdim.begin(), dims.begin(),
                dims.begin() + nc);  // first nc dims

    NC res(ndim);  // ndim.size == sdim.size

    // assign slice of combs applying fc func.
    int cr = 0;  // res comb counter
    auto combRes = combinations(res.dims);
    for (auto &c : combinations(ndim))
      res.assign(fc(slice(c), res.dims), combRes[cr++ % combRes.size()]);

    return res;
  }

 private:
  vector<T> __fromString(const string &s, const string regExp) {
    regex re(regExp);

    vector<T> v;
    for (auto it = sregex_token_iterator(s.begin(), s.end(), re);
         it != sregex_token_iterator(); it++) {
      try {
        v.push_back(stod(*it));
      } catch (const std::invalid_argument &e) {
        cerr << *it << ":invalid number" << endl;
      }
    }

    return v;
  }

#define __regExp R"([-+]?\d*\.?\d+([eE][-+]?\d+)?)"

 public:
  NC fromString(const string &s, const string regExp = __regExp) {
    data = __fromString(s, regExp);

    dims = {(int)data.size()};  // 1d size
    __recalc();

    return *this;
  }

  NC fromFile(const string fname, const string regExp = __regExp) {
    ifstream file(fname);

    if (file.is_open()) {
      data.clear();

      for (string line; getline(file, line);) {
        auto v = __fromString(line, regExp);
        data.insert(data.end(), v.begin(), v.end());
      }

      dims = {(int)data.size()};  // 1d size
      __recalc();
    }
    return *this;
  }

 private:
  // misc
  string fmt(const char *format, ...) {
    string buffer(1024, 0);
    va_list args;
    va_start(args, format);
    vsnprintf((char *)buffer.c_str(), buffer.size(), format, args);
    va_end(args);
    return buffer;
  }

  VVInt combinations(const VInt &inArr) {
    VInt currComb, wArr;

    function<VVInt(int)> genCombs = [&](int ix) -> VVInt {
      VVInt res;

      if (ix == wArr.size())
        res = {currComb};
      else {
        for (int i = 0; i < wArr[ix]; i++) {
          currComb[ix] = i;
          for (auto &c : genCombs(ix + 1)) res.push_back(c);
        }
      }
      return res;
    };

    wArr = inArr.empty() ? dims : inArr;
    currComb.resize(wArr.size());

    return genCombs(0);
  }

 public:
  string dimsStr() {
    string s;
    for (int i = 0; i < ndims; i++) {
      s += to_string(dims[i]);
      if (i < ndims - 1) s += ',';
    }
    return s;
  }
  VInt descToDims(string desc) {  // (#,#,#,#...)
    auto start = desc.find('('), end = desc.find(')', start);
    desc = desc.substr(start + 1, end - start - 1);

    auto _desc = desc.c_str();

    return strToDims(desc);
  }
  VInt strToDims(string desc) {  // #,#,...,#
    size_t start = 0, end = desc.find(',');
    VInt vdim;
    auto _desc = desc.c_str();

    while (end != string::npos) {
      vdim.push_back(stoi(desc.substr(start, end - start)));
      start = end + 1;
      end = desc.find(',', start);
    }
    vdim.push_back(stoi(desc.substr(start)));
    return vdim;
  }

  // numpy interface
  void save(string name) {
    FILE *fh;

    const char *nm = name.c_str();
    if ((fh = fopen(nm, "wb")) != nullptr) {
      string desc =
          fmt("{'descr':'<f%d', 'fortran_order':False, 'shape':(%s,)}",
              sizeof(T), dimsStr().c_str());

      int szHdr = sizeof(HeaderNPY) + desc.length() + 1;  // header size + $0a
      desc += string(64 * ((szHdr / 64) + 1) - szHdr, ' ') + '\x0a';

      HeaderNPY hdr(desc.size());

      // write header, desc, data
      fwrite(&hdr, sizeof(HeaderNPY), 1, fh);
      fwrite(desc.data(), 1, desc.length(), fh);
      fwrite(data.data(), 1, sizeBytes, fh);

      fclose(fh);
    }
  }

  NC load(string name) {
    FILE *fh;
    NC nc;

    if ((fh = fopen(name.c_str(), "rb")) != nullptr) {
      HeaderNPY hdr;
      assert(fread(&hdr, sizeof(hdr), 1, fh) == sizeof(hdr));

      string desc(hdr.headLen, ' ');
      assert(fread(desc.data(), 1, desc.size() == desc.size(), fh));

      nc = NC(descToDims(desc));  // (#,#,#)
      assert(fread(nc.data.data(), 1, nc.sizeBytes, fh) == nc.sizeBytes);

      *this = nc;

      fclose(fh);
    }
    return nc;
  }

  // lin alg
 private:
  void assertQuadratic() {
    assert((ndims >= 2 && dim0 == dim(1)) && "non quadratic matrix");
  }
  void assert_nxn() {
    assert((ndims == 2 && dim0 == dim(1)) && "non N x N matrix");
  }

  inline int &dim(int ix) {  // dims in reverse order
    assert(ix < ndims && "dims index out of range");
    return dims[ndims - 1 - ix];
  }

 public:
  NC coFactor(int x) {  // cofactor remove row 1, current col(x)
    int n = dim0;
    NC res(n - 1, n - 1);

    for (int i = 1, si = 0; i < n; i++, si++)  // remove row(1)
      for (int j = 0, sj = 0; j < n; j++)
        if (j != x) res.at(si, sj++) = at(i, j);  // remove col(x)

    return res;
  }

  T detnxnMT() {
    assertQuadratic();
    return _detnxnMT(*this);
  }

  T _detnxnMT(NC a) {  // cofactor based, very slow !!
    int n = a.dim0;

    switch (n) {
      case 1:
        return a.at(0);
      case 2:
        return a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);  //  Sarrus rule
      case 3:
        return (a(0, 0) * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))) -
               (a(0, 1) * (a(1, 0) * a(2, 2) - a(1, 2) * a(2, 0))) +
               (a(0, 2) * (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0)));
      default: {
        T d = 0;
        for (int x = 0; x < n; x++)
          d += (x % 2 == 0 ? 1 : -1) * a.at(0, x) * a._detnxnMT(a.coFactor(x));

        return d;
      }
    }
  }

  T detBareiss() {
    T res = (T)(1 + 1e-10);
    NC a(*this);

    assertQuadratic();

    int n = dim0;

    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        if (a(i, i) != static_cast<T>(0)) {
          T factor = a(j, i) / a(i, i);
          for (int k = i; k < n; k++) a(j, k) -= factor * a(i, k);
        }
      }
      res *= a(i, i);
    }
    return res;
  }

  NC slice(const VInt &index) {
    assert(!index.empty() && "slice with empty array not allowed");

    auto widx = index;
    widx.insert(widx.end(), ndims - widx.size(), 0);
    int istart = calcIndex(widx);

    widx = index;
    for (int i = index.size(); i < ndims; i++) widx.push_back(dims[i] - 1);
    auto iend = calcIndex(widx) + 1;

    auto sz = ndims - index.size();
    copy(dims.begin() + index.size(), dims.begin() + index.size() + sz,
         widx.begin());  // widx=dims[szIndex..szDims-szIndex]
    widx.resize(sz);

    NC res(widx);

    copy(data.begin() + istart, data.begin() + istart + (iend - istart),
         res.data.begin());

    return res;
  }

  NC det() {  // N x N, MT version
    NC res;

    assertQuadratic();

    if (ndims == 2) {
      res = NC(1);
      res.at(0) = detLU_MT();
    } else {
      VInt tdim(dims.begin() + 0, dims.end() - 2);  // dims[0..-2]
      VVInt tdims = combinations(tdim);
      res = NC(tdim);

      MyThread().runSeq([&](int ix, int nth) {
        for (; ix < tdims.size(); ix += nth)
          res.at(tdims[ix]) = slice(tdims[ix]).detLU_MT();
      });
    }
    return res;
  }

  NC detST() {  // ST version
    NC res;

    assertQuadratic();

    if (ndims == 2) {
      res = NC(1);
      res.at(0) = detBareiss();
    } else {
      VInt tdim(dims.begin() + 0, dims.end() - 2);  // dims[0..-2]
      res = NC(tdim);

      for (auto &t : combinations(tdim)) res.at(t) = slice(t).detBareiss();
    }
    return res;
  }

  void rangeSlice(const VInt &index, int &aStart, int &aEnd) {
    if (index.empty()) {
      aStart = aEnd = 0;
    } else {
      auto widx = index;  // start
      widx.insert(widx.end(), ndims - widx.size(), 0);
      aStart = calcIndex(widx);

      widx = index;  // end
      for (int i = index.size(); i < ndims; i++) widx.push_back(dims[i] - 1);
      aEnd = calcIndex(widx) - 1;
    }
  }

  void assign(const NC a, const VInt &_dims) {
    int aStart, aEnd;
    rangeSlice(_dims, aStart, aEnd);
    for (int i = 0; i < a.size; i++) data[i + aStart] = a.data[i];
  }
  void swap(int i, int j, int k, int l) {
    auto tmp = at(i, j);
    at(i, j) = at(k, l);
    at(k, l) = tmp;
  }
  int prod(VInt &v) {
    int p = 1;
    for (auto i : v) p *= i;
    return p;
  }

  T diagProd() {  // diag. prod. N x N matrix
    T res = (T)1;
    for (int i = 0; i < dim0; i++) res *= at(i, i);
    return res;
  }

  VInt revVInt(VInt a) {  // reserve array
    reverse(a.begin(), a.end());
    return a;
  }

  int transposedCoord(int index) {
    int res = 0;
    for (int i = 0; i < ndims && index != 0; i++) {
      res += (index % dims[i]) * mlt[i];
      index /= dims[i];
    }
    return res;
  }
  NC reshape(VInt ndim) {
    assert(size = prod(ndim) && "reshape: incompatible sizes");
    NC res(ndim);

    res.data = data;
    return res;
  }

  int _min(int a, int b) { return a < b ? a : b; }

  NC dotST(NC &a, NC &b) { return a.dotST(b); }
  NC dot(NC &a, NC &b) { return a.dot(b); }

 private:
  NC __prepDot(NC &a, int &pivotPos, int &pivd, int &prs) {
    pivotPos = _min(1, a.ndims - 1);  // position in a of pivot dim in  1|0
    pivd = a.dim(pivotPos);           // pivot dimension, a.dim(1|0) = dim(0)

    assert(dim0 == pivd && "dot: shapes not aligned");

    VInt sDim(dims.begin(), dims.end() - 1);  // all but dim(0)
    auto aDim = a.dims;                       // remove 'pp' 0|1
    aDim.erase(aDim.begin() + a.ndims - 1 - pivotPos);

    VInt resDim = sDim;  // resDim = sDim + aDim
    resDim.insert(resDim.end(), aDim.begin(), aDim.end());

    if (resDim.empty()) resDim = {1};

    prs = prod(sDim);

    return NC(resDim);
  }

 public:
  NC dotST(NC &a) {
    int pivotPos, pivd, prs;
    NC res = __prepDot(a, pivotPos, pivd, prs);

    switch (pivotPos) {
      case 0: {
        for (int ixs = 0, sStart = 0; ixs < prs; ixs++) {
          T p = 0;  // calc dot 1x1
          for (int r = 0; r < pivd; r++) p += at(sStart + r) * a.at(r);
          res.at(ixs) = p;
          sStart += pivd;  // next self pivd slice position
        }
      } break;

      case 1: {
        VInt raDim;
        if (a.ndims > 2)
          raDim.insert(raDim.begin(), a.dims.begin(), a.dims.end() - 2);

        int ahi1 = (a.ndims > 1) ? a.dim(1) : 1;
        int pra = prod(raDim), aStride = a.size / pra, adim0 = a.dim0,
            ahi0 = a.dim0;

        for (int ixs = 0, ixr = 0, sStart = 0; ixs < prs;
             ixs++, sStart += pivd) {  // 1x2 reduction of self
          for (int ixa = 0, aStart = 0; ixa < pra;
               ixa++, aStart += aStride) {  // next a slice

            for (int c = 0; c < ahi0; c++) {
              T p = 0;
              for (int r = 0; r < ahi1; r++)
                p += at(sStart + r) * a.at(aStart + r * adim0 + c);
              res.at(ixr++) = p;
            }
          }
        }
      } break;

      default:
        assert("dot: internal error pivot not 0 or 1");
    }

    return res;
  }

  NC dot(NC &a) {
    int pivotPos, pivd, prs;
    NC res = __prepDot(a, pivotPos, pivd, prs);

    switch (pivotPos) {
      case 0: {
        MyThread().runSeq([&](int ix, int nth) {
          for (int sStart = ix * pivd; ix < prs; ix += nth) {
            T p = 0;
            for (int r = 0; r < pivd; r++) p += at(sStart + r) * a.at(r);
            res.at(ix) = p;
          }
        });
      } break;
      case 1: {
        VInt raDim;
        if (a.ndims > 2)
          raDim.insert(raDim.begin(), a.dims.begin(), a.dims.end() - 2);

        int pra = prod(raDim), aStride = a.size / pra, adim0 = a.dim0,
            ahi0 = a.dim0;
        int ahi1 = (a.ndims > 1) ? a.dim(1) : 1;

        MyThread().runSeq([&](int ix, int nth) {
          for (; ix < prs; ix += nth) {
            for (int ixa = 0, aStart = 0, sStart = ix * pivd,
                     ixr = ix * (pra * adim0);
                 ixa < pra; ixa++, aStart += aStride) {
              for (int c = 0; c < ahi0; c++) {
                T p = 0;
                for (int r = 0; r < ahi1; r++)
                  p += at(sStart + r) * a.at(aStart + r * adim0 + c);
                res.at(ixr++) = p;
              }
            }
          }
        });
      } break;
      default:
        assert("dot: internal error pivot not 0 or 1");
    }

    return res;
  }

  NC inv_nxn() {
    int n = dim0;
    NC res = NC::identity(n);

    NC a(*this);

    for (int j = 0; j < n; j++) {
      for (int i = j; i < n; i++) {
        if (a(i, j) != 0) {
          for (int k = 0; k < n; k++) {
            a.swap(j, k, i, k);
            res.swap(j, k, i, k);
          }
          auto tmp = 1 / a(j, j);
          for (int k = 0; k < n; k++) {
            a(j, k) = tmp * a(j, k);
            res(j, k) = tmp * res(j, k);
          }

          for (int k = 0; k < n; k++) {
            if (k != j) {
              tmp = -a(k, j);
              for (int c = 0; c < n; c++) {
                a(k, c) = a(k, c) + tmp * a(j, c);
                res(k, c) = res(k, c) + tmp * res(j, c);
              }
            }
          }
          break;
        }
      }
    }

    return res;
  }

  NC invST() {  // matrix inversion ST
    assertQuadratic();

    NC res(*this);
    if (ndims == 2)
      res = inv_nxn();
    else {
      VInt tdim(dims.begin() + 0, dims.end() - 2);  // dims[0..-2]
      for (auto &t : combinations(tdim)) res.assign(slice(t).inv_nxn(), t);
    }
    return res;
  }

  NC inv() {  // MT version
    NC res;

    assertQuadratic();

    if (ndims == 2) {
      res = NC(1);
      res.at(0) = detBareiss();
    } else {
      VInt tdim(dims.begin() + 0, dims.end() - 2);  // dims[0..-2]
      VVInt tdims = combinations(tdim);
      res = *this;

      MyThread().runSeq([&](int ix, int nth) {
        for (; ix < tdims.size(); ix += nth)
          res.assign(slice(tdims[ix]).inv_nxn(), tdims[ix]);
      });
    }
    return res;
  }

  NC solve(NC &b) { return inv().dot(b); }
  NC solveST(NC &b) { return invST().dotST(b); }

  NC transposeST() {
    NC res(revVInt(dims));

    for (int i = 1; i < size / 2 - 1; i++)
      res.at(res.transposedCoord(i)) = at(i);
  }

  NC transpose() {  // MT version
    NC res(revVInt(dims));

    MyThread(size / 2 - 1).run([&](int t, int from, int to) {
      for (int i = from; i < to; i++) res.at(res.transposedCoord(i)) = at(i);
    });
  }

  T detLU() {  // a bit slower than bareiss
    assert_nxn();

    NC L(dims), U(dims);  // LU decomp.
    int n = dim0;

    for (int i = 0; i < n; i++) {  // L
      for (int j = 0; j < n; j++) {
        if (j >= i) {
          L(j, i) = at(j, i);
          for (int k = 0; k < i; k++) {
            L(j, i) -= L(j, k) * U(k, i);
          }
        }
      }
      for (int j = 0; j < n; j++) {  // U
        if (j >= i) {
          if (j == i)
            U(i, j) = 1;
          else {
            U(i, j) = at(i, j) / L(i, i);
            for (int k = 0; k < i; k++) {
              U(i, j) -= (L(i, k) * U(k, j)) / L(i, i);
            }
          }
        }
      }
    }

    T detL = 1, detU = 1;  // mult L*U
    for (int i = 0; i < n; i++) {
      detL *= L.at(i, i);
      detU *= U.at(i, i);
    }
    return detL * detU;
  }

  T detLU_MT() {  // fastest method

    auto Th = MyThread(dim0);
    barrier sync_point(Th.nth);

    NC a(*this);  // work copy

    Th.run([&](int t, int start, int end) {
      for (int k = 0; k < dim0; ++k) {
        sync_point.arrive_and_wait();  // Synchronize threads

        for (int i = start; i < end; ++i)
          if (i > k) {
            a(i, k) /= a(k, k);
            for (int j = k + 1; j < dim0; ++j) a(i, j) -= a(i, k) * a(k, j);
          }
      }
    });

    return a.diagProd();
  }

  T detLU_01() {
    NC a(*this);  // work copy
    int n = dim0;

    for (int k = 0; k < n; ++k) {
      for (int i = k + 1; i < n; ++i) {
        a(i, k) /= a(k, k);
        for (int j = k + 1; j < n; ++j) a(i, j) -= a(i, k) * a(k, j);
      }
    }

    return a.diagProd();
  }

  T detBareissMT() {  //  n x n
    assert_nxn();

    auto Th = MyThread(dim0);
    barrier sync_point(Th.nth);

    NC a(*this);  // working copy

    Th.run([&](int t, int start, int end) {
      for (int k = 0; k < a.dim0; ++k) {
        sync_point.arrive_and_wait();  // Synchronize threads at the start of
                                       // each step

        for (int i = start; i < end; ++i) {
          if (i <= k)
            continue;  // Skip rows that are not part of the current step

          for (int j = k + 1; j < a.dim0; ++j) {
            // lock_guard<mutex> lock(mtx); // not required, no race conditions
            if (k == 0) {
              a(i, j) = a(i, j) * at(k, k) - at(i, k) * a(k, j);
            } else {
              a(i, j) =
                  (a(i, j) * a(k, k) - a(i, k) * a(k, j)) / a(k - 1, k - 1);
            }
          }
        }
      }
    });
    return a.data.back();  // det is last item a(n-1,n-1)
  }
};
