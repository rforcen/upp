//
//  Timer.h
//  Vect
//
//  Created by asd on 17/03/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#ifndef Timer_h
#define Timer_h

//#include <ctime>
#include <chrono>
#include <functional>
#include <ratio>
#include <thread>

class Timer {
  std::chrono::high_resolution_clock::time_point begin, end;

public:
  Timer() { start(); }
  void start() { begin = std::chrono::high_resolution_clock::now(); }
  long lap() {
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span = end - begin;
    return long(time_span.count());
  }
  long timer(std::string m = "") {
    auto l = lap();
    printf("%s", m.c_str());
    printf(": %ld\n", l);
    start();
    return l;
  }
  long chrono(std::function<void(void)> const &lambda) {
    start();
    lambda();
    return lap();
  }
  static void pause(int secs) {
    std::this_thread::sleep_for(std::chrono::seconds(secs));
  }
};

#endif /* Timer_h */
