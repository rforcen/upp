g++ -shared -static -o mandelbrot.dll -O3 -std=c++2a -I. colormap.c mandel.cpp mandel32.cpp mandel64.cpp mandel128.cpp
copy mandelbrot.dll ..