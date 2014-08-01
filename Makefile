CFLAGS=-std=c++0x -g
LDFLAGS=-lboost_program_options -lboost_system -lboost_filesystem 

OBJECTS = main.o \
  Support.o

all: main 

main: $(OBJECTS)
	g++ $(OBJECTS) -o $@ $(LDFLAGS) 

main.o: main.cpp Support.h Header.h
	g++ -c $(CFLAGS) $< -o $@

Support.o: Support.cpp Support.h Header.h
	g++ -c $(CFLAGS) $< -o $@
clean:
	rm -f *.o *~ main gmon.out 

