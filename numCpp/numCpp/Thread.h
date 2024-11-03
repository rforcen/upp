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
  FromTo(size_t from, size_t to) : from(from), to(to) {}
  size_t from, to;
};

class MyThread {
 public:
  MyThread() : nth(nCpus()), threads(vector<thread>(nth)), size(nCpus()) {}
  MyThread(size_t size)
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
    size_t stride = size / nth, first = 0, extra = size % nth;

    ranges.clear();
    for (size_t i = 0; i < nth; i++) {
      size_t last = first + stride;
      if (extra > 0) {
        extra--;
        last++;
      }
      ranges.push_back(FromTo(first, last));
      first = last;
    }
    // for (auto r : ranges)  printf("%d, ", r.second - r.first);
  }

  size_t from(int t) { return ranges[t].from; }
  size_t to(int t) { return ranges[t].to; }

  void run(std::function<void(int, size_t, size_t)> const &lambda) {  // t, from, to
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() { lambda(t, from(t), to(t)); });
    }
    for (int t = 0; t < nth; t++) threads[t].join();
  }

  void run(std::function<void(size_t)> const &lambda) {  // i
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() {
        for (size_t i = from(t); i < to(t); i++) lambda(i);
      });
    }
    for (int t = 0; t < nth; t++) threads[t].join();
  }

  void run(std::function<void(int, size_t)> const &lambda) {  // t, i
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
        for (size_t i = from(t); i < to(t); i++) lambda();
      });
    }
    for (int t = 0; t < nth; t++) threads[t].join();
  }

  int nth = nCpus();
  size_t size = 0;
  vector<thread> threads;
  vector<FromTo> ranges;
};

#endif /* Thread_h */

