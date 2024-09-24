all:
	@g++ -c -O3 mandel32.cpp mandel64.cpp mandel128.cpp
	@g++ -c mandelmpfr.cpp
	@ar -r libmandelbrot.a *o
	@mv libmandelbrot.a ..
	@rm -fr *.o
