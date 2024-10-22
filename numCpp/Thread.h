//
//  MyThread.h
//

#ifndef Thread_h
#define Thread_h

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <utility>
#include <functional>
#include <thread>

using std::thread;
using std::vector;

class FromTo {
 public:
  FromTo(int from, int to) : from(from), to(to) {}
  int from, to;
};

class MyThread {
 public:
  MyThread() : nth(nCpus()), threads(vector<thread>(nth)), size(nCpus()) {}
  MyThread(int size)
      : nth(size < nCpus() ? size : nCpus()),
        size(size),
        threads(vector<thread>(nth)) {
    chunk_ranges();
  }

  void runSeq(std::function<void(int, int)> const &lambda) {  // seq
    for (int t = 0; t < nth; t++)
      threads[t] = thread([this, lambda, t]() { lambda(t, nth); });

    for (int t = 0; t < nth; t++) threads[t].join();
  }

  static int nCpus() { return int(thread::hardware_concurrency()); }

  void chunk_ranges()  // Create the result and calculate the stride size and
                       // the remainder if any.
  {
    int stride = size / nth, first = 0, extra = size % nth;

    ranges.clear();
    for (int i = 0; i < nth; i++) {
      int last = first + stride;
      if (extra > 0) {
        extra--;
        last++;
      }
      ranges.push_back(FromTo(first, last));
      first = last;
    }
    // for (auto r : ranges)  printf("%d, ", r.second - r.first);
  }

  int from(int t) { return ranges[t].from; }
  int to(int t) { return ranges[t].to; }

  void run(std::function<void(int, int, int)> const &lambda) {  // t, from, to
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() { lambda(t, from(t), to(t)); });
    }
    for (int t = 0; t < nth; t++) threads[t].join();
  }

  void run(std::function<void(int)> const &lambda) {  // i
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() {
        for (int i = from(t); i < to(t); i++) lambda(i);
      });
    }
    for (int t = 0; t < nth; t++) threads[t].join();
  }

  void run(std::function<void(int, int)> const &lambda) {  // t, i
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() {
        for (int i = from(t); i < to(t); i++) lambda(t, i);
      });
    }
    for (int t = 0; t < nth; t++) threads[t].join();
  }

  void run(std::function<void(void)> const &lambda) {  // ()
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() {
        for (int i = from(t); i < to(t); i++) lambda();
      });
    }
    for (int t = 0; t < nth; t++) threads[t].join();
  }

  int nth = nCpus(), size = 0;
  vector<thread> threads;
  vector<FromTo> ranges;
};

#endif /* Thread_h */

