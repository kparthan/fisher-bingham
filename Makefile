# kernel-compatible version
#CFLAGS=-std=c++98 -c -I/home/parthan/external_libs/ -fopenmp
#LDFLAGS=-static -lboost_program_options -lboost_filesystem -fopenmp -lnlopt -lm

CFLAGS=-std=c++0x -c -O3 -fopenmp
#CFLAGS=-std=c++0x -g -c -fopenmp
LDFLAGS=-lboost_program_options -lboost_system -lboost_filesystem -fopenmp -lnlopt -lm

OBJECTS = main.o \
  Support.o \
  Normal.o  \
  vMC.o \
  vMF.o \
  FB4.o \
  Kent.o \
  FB6.o \
  Optimize.o \
  Optimize2.o \
  Structure.o \
  Mixture.o \
  Test.o \
  Experiments.o

all: main 

main: $(OBJECTS)
	g++ $(OBJECTS) -o $@ $(LDFLAGS) 

main.o: main.cpp 
	g++ $(CFLAGS) $< -o $@

Support.o: Support.cpp Support.h Header.h UniformRandomNumberGenerator.h
	g++ $(CFLAGS) $< -o $@

Normal.o: Normal.cpp Normal.h 
	g++ $(CFLAGS) $< -o $@

vMF.o: vMF.cpp vMF.h Header.h
	g++ $(CFLAGS) $< -o $@

vMC.o: vMC.cpp vMC.h Header.h
	g++ $(CFLAGS) $< -o $@

FB4.o: FB4.cpp FB4.h Header.h
	g++ $(CFLAGS) $< -o $@

Kent.o: Kent.cpp Kent.h Header.h
	g++ $(CFLAGS) $< -o $@

FB6.o: FB6.cpp FB6.h Header.h
	g++ $(CFLAGS) $< -o $@

Structure.o: Structure.cpp Structure.h Header.h
	g++ $(CFLAGS) $< -o $@

Mixture.o: Mixture.cpp Mixture.h Header.h
	g++ $(CFLAGS) $< -o $@

Optimize.o: Optimize.cpp Optimize.h Header.h
	g++ $(CFLAGS) $< -o $@

Optimize2.o: Optimize2.cpp Optimize2.h Header.h
	g++ $(CFLAGS) $< -o $@

Test.o: Test.cpp Test.h Header.h
	g++ $(CFLAGS) $< -o $@

Experiments.o: Experiments.cpp Experiments.h Header.h
	g++ $(CFLAGS) $< -o $@

clean:
	rm -f *.o *~ main gmon.out 

