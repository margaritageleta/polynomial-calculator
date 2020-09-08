CXXFLAGS = -Wall -std=c++11 -O2

all: main.exe

clean:
	rm -f main.exe *.o

main.exe: main.o Polynomial.o
	$(CXX) $^ -o $@

main.cc: Polynomial.hh

Polynomial.cc: Polynomial.hh