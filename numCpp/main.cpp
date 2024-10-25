/*
 NC test
*/
#include <iostream>
#include <math.h>

#include "Timer.h"
#include "Thread.h"
#include "numCpp.h"

// #define _USE_BOOST
// #define _USE_MPREAL
// #define _USE_HALF

typedef long double Real;

void test01() {
  puts("test01");

  NC<double> n1;
  NC<float> n2(3, 4, 5);
  NC<float> n21;
  auto n3 = NC<float>::random(3, 4, 6);
  auto n4 = NC<float>::random(4, 5, 6);

  auto n5 = NC<float>::arange(100);

  n21 = n3 - n3;

  n4.print("n4");
  n5.print("n5");

  n21.print("n21");

  auto v = n3.at(2, 3, 5);
  auto v1 = n3.at(2, 3, 5);

  n3.at(2, 2, 2) = v;
  if (n3.at(2, 2, 2) == v)
    cout << "index ok" << endl;
  else
    cout << "index wrong" << endl;

  cout << v << endl << v1 << endl;

  cout << endl << "sum:" << n3.sum() << endl;
  cout << endl << "sum:" << n4.sum() << endl;
  cout << endl << "sum:" << n21.sum() << endl;

  cout << "--------- n3 ------------" << endl;

  n3.print("n3");
  cout << endl;
  n3.save("n3");

  cout << "--------- n4 ------------" << endl;
  n4.load("n3");
  n4.print("n4");
  cout << endl;
  cout << (!(n3 != n4) ? "n3==n4" : "n3!=n4") << endl;
  cout << ((n3 == n4) ? "n3==n4" : "n3!=n4") << endl;
}

void test02() {
  /*
NC<double> n1;
for (auto &c : n1.combinations(4, 2, 3)) {
for (auto &i : c) cout << i << ' ';
cout << endl;
}*/

  puts("test02");

  auto n2 = NC<double>::random(400, 400);
  n2.save("n2");
  cout << endl << n2.detBareiss() << endl;
}

void test03() {
  puts("test03");

  auto n0 = NC<float>::random(3, 4, 6);
  auto n1 = n0.slice(1, 1);
  n0.print();
  cout << endl;
  n1.print();
  cout << endl;
}

void test04() {
  puts("test04");

  auto n0 = NC<double>::random(30, 5, 400, 400);

  NC d;

  for (int i = 0; i < 3; i++) {
    n0.rand();
    auto lap = Timer().chrono([&] { d = n0.det(); });
    cout << "lap:" << lap << "ms" << endl;
  }

  /*n0.print();
  cout << endl;*/

  d.print();
  cout << endl;

  n0.save("n0");
  d.save("d");
}

void test05() {
  puts("test05");

  auto n0 = NC<double>::random(3, 5, 40, 40);
  auto d = n0.detST();

  n0.save("n0");
  d.save("d");
}

void test06() {
  puts("test06");
  auto n0 = NC<double>::random(3, 6, 40, 40);
  auto n0Inv = n0.inv();

  n0.save("n0");
  n0Inv.save("n0Inv");
}

void test07() {
  puts("test07");
  int pivot = 15;
  auto n0 = NC<double>::random(90, 90, pivot);
  auto n1 = NC<double>::random(90, pivot, 90);
  NC nDot;

  for (int i = 0; i < 3; i++) {
    cout << "lapMT:" << Timer().chrono([&] { nDot = n0.dot(n1); }) << endl;
    cout << "lapST:" << Timer().chrono([&] { nDot = n0.dotST(n1); }) << endl;
  }

  n0.save("n0");
  n1.save("n1");
  nDot.save("nDot");
}

void test070() {
  puts("test070");
  int pivot = 15;
  auto n0 = NC<double>::random(90, 90, pivot);
  auto n1 = NC<double>::random(90, pivot, 90);
  NC nDot;

  cout << "lapST:" << Timer().chrono([&] { nDot = n0.dotST(n1); }) << endl;
  n0.save("n0");
  n1.save("n1");
  nDot.save("nDot");
}

void test08() {
  puts("test08");
  auto n0 = NC<long double>::random(11, 11);
  long double d = 0, dm = 0;

  cout << "lap detST:" << Timer().chrono([&] { d = n0.detBareiss(); }) << " ms"
       << endl;
  cout << "lap detMT:" << Timer().chrono([&] { dm = n0.detnxnMT(); }) << " ms"
       << endl;

  cout << d << ' ' << dm << ' ' << abs(d - dm) << endl;
}

void test09(size_t n = 400) {
  puts("test09");
  auto n0 = NC<long double>::random(n, n);

  long double d = 0, dm = 0, db = 0, dlumt;

  cout << "lap detMT  :" << Timer().chrono([&] { db = n0.detBareissMT(); })
       << " ms,\t\tdet:" << db << endl;
  cout << "lap detST  :" << Timer().chrono([&] { d = n0.detBareiss(); })
       << " ms,\t\tdet:" << d << endl;
  cout << "lap detLU  :" << Timer().chrono([&] { dm = n0.detLU(); })
       << " ms,\t\tdet:" << dm << endl;
  cout << "lap detLUMT:" << Timer().chrono([&] { dlumt = n0.detLU_MT(); })
       << " ms,\t\tdet:" << dlumt << endl;

  // cout << d << ' ' << dm << ' ' << db << endl;
  cout << "relative diff:" << abs(d - dm) / d << endl;
}

#ifdef _USE_BOOST
// boost decl. best known option for MP
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/format.hpp>

using namespace boost::multiprecision;
using namespace boost;
// ----------

void test10(size_t n = 400) {
  puts("test10");
  using fpb = number<cpp_bin_float<20>>;  // 64 bits = 20 digits, 32
                                          // bit exponent  each digit
                                          // is about 3.32 bits

  // printf format a cpp_dec
  auto bfmt = [=](fpb f, const char* fmt) { return str(format(fmt) % f); };

  auto n1 = NC<long double>::random(n, n);
  NC<fpb> n2(n1);

  long double dn1 = 0;
  fpb dn2 = 0;

  printf("det's for n:%d\n", n);

  auto l1 = Timer().chrono([&] {
    dn1 = n1.detLU_MT();
  });  // boost mp is about 33 times slower than c++ long double
  auto l2 = Timer().chrono([&] { dn2 = n2.detLU_MT(); });

  printf("lap det MT long double      :%5ld ms, det :%8Lg\n", l1, dn1);
  printf("lap det MTboost             :%5ld ms, det :%8s\n", l2,
         bfmt(dn2, "%8Lg").c_str());
  auto relDiff = abs(dn1 - dn2) / dn2;
  printf("rel. diff: %s\n", bfmt(relDiff, "%f").c_str());
}
#endif

#ifdef _USING_MPREAL  // slower than boost
#include "mpreal.h"
using mpfr::mpreal;

void test11(size_t n = 300) {
  puts("test11");
  mpreal::set_default_prec(mpfr::digits2bits(20));

  auto n1 = NC<long double>::random(n, n);
  NC<mpreal> n2(n1);

  long double dn1 = 0;
  mpreal dn2 = 0;

  printf("det's for n:%d\n", n);

  auto l1 = Timer().chrono([&] {
    dn1 = n1.detLU_MT();
  });  // boost mp is about 33 times slower than c++ long double
  auto l2 = Timer().chrono([&] { dn2 = n2.detLU_MT(); });

  printf("lap det MT long double      :%5ld ms, det :%8Lg\n", l1, dn1);
  printf("lap det MT mpreal           :%5ld ms, det :%8Lg\n", l2,
         dn2.toLDouble());
  auto relDiff = abs(dn1 - dn2) / dn2;
  printf("rel. diff: %s\n", relDiff.toString("%Lg").c_str());
}
#endif

void test12(int n = 400) {
  puts("test12");
  auto n1 = NC<long double>::random(n, n, n);
  auto l1 = Timer().chrono([&] { n1 = n1.sort(); });

  n1.print();
  printf("matrix size:%ld, lap sort :%5ld ms\n", n1.size, l1);
}

void test13(size_t n = 400) {
  puts("test13");
  auto n1 = NC<long double>::random(n);
  n1.print();
  auto l1 = Timer().chrono([&] { n1 += 1; });
  n1.print();
  auto n2 = n1 + n1;
  n2.print();

  auto n3 = (n1 + n2) + n1;
  n3.print();
  n3 /= 2;
  n3.print();

  n1.inApply([=](long double x) { return sin(x); });
  n1.print();

  printf("matrix size:%ld, lap sort :%5ld ms\n", n1.size, l1);
}

void test14(int n = 400) {
  puts("test14");
  auto n1 = NC<long double>::random(n, n, n);
  auto n2(n1);

  auto l1 = Timer().chrono(
      [&] { n1.inApply([=](long double x) { return sin(x); }); });

  printf("matrix size:%ld\nlap apply   :%5ld ms\n", n1.size, l1);
  l1 = Timer().chrono(
      [&] { n2.inApplyMT([=](long double x) { return sin(x); }); });

  printf("lap apply MT:%5ld ms\nresult:%d\n", l1, n1 == n2);
}

void test15(size_t n = 400) {  // << oper -> append
  puts("test15");
  NC n1 = NC<double>::random(n, n);
  auto n2(n1);

  printf(
      "n1==n2:%d, n1!=n2:%d, n1>=n2:%d, n1>n2:%d, n1<n2:%d, n1<=n2:%d, "
      "\n\n",
      n1 == n2, n1 != n2, n1 >= n2, n1 > n2, n1 < n2, n1 <= n2);

  n1 << n2;

  n1.print("n1");
  n2.print("n2");

  for (int i = 0; i < n; i++) n1 << NC<double>::random(n);

  n1.print("n1");

  auto c = n1.concat(n1).reshape({n + 2, n, n});
  c.print("c");
}

void test16() {  // push_back aka <<
  puts("test16");
  using nfp = NC<double>;

  (nfp(4).lin() << nfp(2).lin(5)).print("a");
  (nfp(4).lin() << nfp(4).lin(5)).print("a");
  (nfp(3, 4).lin() << nfp(4).lin(15)).print("a");
  (nfp(2, 3, 4).lin() << nfp(2, 3, 4).lin(35)).print("a");
  (nfp(2, 3, 4).lin() << nfp(3, 4).lin(35)).print("a");

  {
    try {
      (nfp(2, 3, 4).lin() << nfp(5, 4).lin(35)).print("a");

    } catch (runtime_error) {
      // printf("\n%s + %s: incompatible dims!!\n", a.dimsStr().c_str(),
      // b.dimsStr().c_str());
      puts("\nincompatible dims!!");
    }
  }

  // case 2 1 dim eq dims
}

void test17(size_t n = 5) {  // filter
  puts("test17");
  NC n1 = NC<double>::random(n, n) + 0.1;

  auto nflt = n1.filter([=](double x) { return x > 0.7; });

  n1.print("n1");
  nflt.print("nflt");

  auto nfs = n1.filter("x > 0.7");
  nfs.print("nflt");
  cout << "eq?:" << (nflt == nfs) << endl;

  n1 = NC<double>::random(10000000);

  auto lapf = Timer().chrono([&] {
    nfs = n1.filter([=](double x) { return (x >= 0.7) & (x <= 0.9); });
  });
  printf("lap filter expr c++  for %ld items:%ld ms\n", n1.size, lapf);

  auto nfcpp = nfs;

  auto lapsf =
      Timer().chrono([&] { nfs = n1.filter("(x >= 0.7) & (x <= 0.9)"); });
  printf("lap filter expr str  for %ld items:%ld ms\n", n1.size, lapsf);

  if (nfcpp == nfs)
    cout << "string eval. OK!" << endl;
  else
    cout << "error in string eval." << endl;
}

void test18(size_t n = 400) {
  puts("test18");
  NC n1 = NC<double>::random(n, n);

  n1 += 10;
  n1.print("n1+10");
  n1 -= 10;
  n1.print("n1-10");
  n1 *= 10;
  n1.print("n1*10");
  n1 /= 10;
  n1.print("n1/10");

  auto m = n1.match([=](double x) { return x > 0.6 && x < 0.8; }).sort();
  m.print("m:");
}

void test19() {
  puts("test19");
  NC n1;

  n1.fromString("1,2,3,4,5 6\t7;8; ;9  10");
  n1.print("n1:");

  string s;
  for (int i = 0; i < 20; i++) s += to_string(i) + " ,;\t"[rand() % 4];
  cout << "fromString:" << s << endl;
  n1.fromString(s);
  n1.print("n1:");

  n1.fromFile("/home/asd/MyApps/numCpp/main.cpp");
  n1.print("n1 frm main.cpp:");
}

void test20() {
  puts("test20");
  NC n1(2, 2);
  n1.lin();

  n1.print("n1:");

  NC nx = n1.expand({2, 3});
  nx.print("nx:");

  auto n2 = nx.contract(2, [=](NC<double> nc, VInt v) {
    NC res(v);

    auto s = nc.sum();
    res.fill(s);

    return res;
  });
  n2.print("n2:");
}

#ifdef _USE_HALF

#include "half.h"
void test21() {
  puts("test21");
  size_t n = 30;
  using half_float::half;
  NC<half> a = NC<half>::random(n, n), b(n, n);
  b.rand();

  cout << "size of half:" << sizeof(half) << endl;  // 2 bytes

  a.slice(n / 3).range(0, n / 3).print("a:");
  b.slice(n / 3).range(0, n / 3).print("b:");

  auto c = a.slice(n / 3, 0);  // at(n/3,0)
  cout << "size of c:" << c.size << endl;
  c.print("c:");

  // compare performance w/float: float:half performance ratio -> 6:70
  // 1:11

  auto lh = Timer().chrono([&] {
    for (int i = 0; i < n * 100; i++) a.detST();
  });
  printf("%d items, lap half:%ld\n", n, lh);

  NC<float> af(a);
  auto lf = Timer().chrono([&] {
    for (int i = 0; i < n * 100; i++) af.detST();
  });

  printf("%d items, lap float:%ld\n", n, lf);
}
#endif

void test22() {
  puts("test22");

  for (int i = 0; i < 8; i++) {
    NC a = NC<double>::random(6, 30, 30, 30);

    // enhace compression ratio by reducing # of different #'s
    a.inApplyMT([&](double x) { return floor(x * 100) / 100.0; });

    auto sa = a.sum();

    Timer t;

    a.compress();

    printf(
        "%2d: org. size:%8ld bytes, compressed size:%8ld bytes, ratio:%.1f %%, "
        "lap:%4ld ms, ",
        i, a.sizeBytes, a.compSize(),
        100 - floor(100.0 * a.compSize() / a.sizeBytes), t.lap());

    a.decompress();

    puts((sa == a.sum()) ? "decompress ok" : "decompression failure!");
  }
}

void test23() {
  auto a = NC<double>::random(10, 10), b(a);

  cout << "a==b : " << (a == b) << endl;

  a.det().print("det a:");
  b.det().print("det b:");
}

void testAll() {
  test01();
  test02();
  test03();

  test04();
  test05();
  test06();
  test07();
  test070();
  test08();
  test09();
  test12();
  test13();
  test14();
  test15();
  test16();
  test17();
  test18();
  test19();
  test20();
  test22();
}

int main(int argc, const char* argv[]) {
  srand(time(nullptr));

  test17();

  return 0;
}

