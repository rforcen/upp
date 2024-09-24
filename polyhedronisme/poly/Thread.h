//
//  Thread.h
//  HillClimbImage
//
//  Created by asd on 29/04/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#ifndef Thread_h
#define Thread_h

#include "common.hpp"
#include <mutex>
#include <thread>

using std::thread, std::mutex;

class Thread {
public:
  Thread()
      : nth(int(thread::hardware_concurrency())), threads(new thread[nth]) {}
  Thread(int size)
      : nth(size < int(thread::hardware_concurrency())
                ? size
                : int(thread::hardware_concurrency())),
        segSz(size > nth ? size / nth : 1), size(size),
        threads(new thread[nth]), mtx(new mutex) {}

  ~Thread() {
    delete[] threads;
    delete mtx;
  }
  static int getnthreads() { return int(thread::hardware_concurrency()); }

  int from(int t) { return t * segSz; }
  int to(int t) { return ((t == nth - 1) ? size : (t + 1) * segSz); }

  void run(std::function<void(int, int, int)> const &lambda) { // t, from, to
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() { lambda(t, from(t), to(t)); });
    }
    for (int t = 0; t < nth; t++)
      threads[t].join();
  }

  void run_once_per_thread(std::function<void(int)> const &lambda) { // i
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([lambda, t]() { lambda(t); });
    }
    for (int t = 0; t < nth; t++)
      threads[t].join();
  }

  void run(std::function<void(int)> const &lambda) { // i
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() {
        for (int i = from(t); i < to(t); i++)
          lambda(i);
      });
    }
    for (int t = 0; t < nth; t++)
      threads[t].join();
  }

  void run(std::function<void(int, int)> const &lambda) { // t, i
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() {
        for (int i = from(t); i < to(t); i++)
          lambda(t, i);
      });
    }
    for (int t = 0; t < nth; t++)
      threads[t].join();
  }

  void run(std::function<void(void)> const &lambda) { // ()
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() {
        for (int i = from(t); i < to(t); i++)
          lambda();
      });
    }
    for (int t = 0; t < nth; t++)
      threads[t].join();
  }

  void run(std::function<void(int, mutex *mtx)> const &lambda) { // i
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() {
        for (int i = from(t); i < to(t); i++)
          lambda(i, mtx);
      });
    }
    for (int t = 0; t < nth; t++)
      threads[t].join();
  }

  void run(std::function<void(int, int, mutex*mtx)> const &lambda) { // t, i
    for (int t = 0; t < nth; t++) {
      threads[t] = thread([this, lambda, t]() {
        for (int i = from(t); i < to(t); i++)
          lambda(t, i, mtx);
      });
    }
    for (int t = 0; t < nth; t++)
      threads[t].join();
  }


  void lock() { mtx->lock(); }
  void unlock() { mtx->unlock(); }

  int nth = int(thread::hardware_concurrency()), segSz = 0, size = 0;
  thread *threads = nullptr;

  mutex *mtx = nullptr; // same mutex for all threads
};

#endif /* Thread_h */
