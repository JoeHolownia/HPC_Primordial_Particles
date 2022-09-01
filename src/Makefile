CFLAGS = -std=c++11 -pg -g -O0 -lm

all: main.o universe.o ioparser.o
	g++ main.o universe.o ioparser.o ${CFLAGS} -o particles

main.o: main.cpp Universe.h IOParser.h
	g++ -c ${CFLAGS} main.cpp

universe.o: Universe.cpp Universe.h Particle.h
	g++ -c ${CFLAGS} Universe.cpp

ioparser.o: IOParser.cpp IOParser.h Particle.h nlohmann/json.hpp
	g++ -c ${CFLAGS} IOParser.cpp
