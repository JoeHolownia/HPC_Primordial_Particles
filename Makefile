CFLAGS = -std=gnu99 -pg -g -O0 -lm

output: main.o universe.o
	g++ main.o universe.o ${CFLAGS} -o decay

main.o: main.cpp Universe.h
	g++ -c ${CFLAGS} main.cpp

universe.o: Universe.cpp Universe.h Particle.h
	g++ -c ${CFLAGS} Universe.cpp
