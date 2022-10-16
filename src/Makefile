MPICXX?=mpic++
CFLAGS = -std=c++11 -O0 -lm

default: all

all: main.o universe.o ioparser.o
	$(MPICXX) main.o universe.o ioparser.o ${CFLAGS} -o particles

main.o: main.cpp universe.h ioparser.h
	$(MPICXX) ${CFLAGS} -c main.cpp

universe.o: universe.cpp universe.h particle.h
	$(MPICXX) ${CFLAGS} -c universe.cpp

ioparser.o: ioparser.cpp ioparser.h particle.h nlohmann/json.hpp
	$(MPICXX) ${CFLAGS} -c ioparser.cpp

clean:
	rm *.o
	rm *.exe
