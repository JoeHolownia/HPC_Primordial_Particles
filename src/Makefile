CFLAGS = -std=c++11 -pg -O0

all: main.o universe.o ioparser.o
	g++ main.o universe.o ioparser.o ${CFLAGS} -o particles

main.o: main.cpp universe.h ioparser.h
	g++ -c ${CFLAGS} main.cpp

universe.o: universe.cpp universe.h particle.h
	g++ -c ${CFLAGS} universe.cpp

ioparser.o: ioparser.cpp ioparser.h particle.h nlohmann/json.hpp
	g++ -c ${CFLAGS} ioparser.cpp

clean:
	rm *.o
	rm *.exe
