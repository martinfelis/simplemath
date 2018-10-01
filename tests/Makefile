CXX=c++
CXXFLAGS=-g

all: main.o catch.hpp tests.cpp ../include/SimpleMath/SimpleMath.h
	${CXX} ${CXXFLAGS} tests.cpp main.o -o runtests

main.o: catch.hpp
	${CXX} -c main.cpp

clean:
	rm -f main.o
	rm -f runtests
