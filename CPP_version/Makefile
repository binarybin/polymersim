polymersim: *.cpp
	g++ -std=c++11 -DNDEBUG -O3 -o $@ *.cpp

no_two: *.cpp
	g++ -std=c++11 -DNDEBUG -DNO_TWO_END -O3 -o polymersim_$@ *.cpp

clean:
	rm polymersim
	rm polymersim_no_two
