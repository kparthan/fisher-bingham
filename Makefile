CFLAGS=-std=c++0x -g 
LDFLAGS=-lboost_program_options -lboost_system -lboost_filesystem -lmpfr

OBJECTS = main.o \
  Support.o \
  Normal.o  \
  vMF.o \
  FB4.o \
  FB6.o \
  Test.o

all: main 

main: $(OBJECTS)
	g++ $(OBJECTS) -o $@ $(LDFLAGS) 

main.o: main.cpp 
	g++ -c $(CFLAGS) $< -o $@

Support.o: Support.cpp Support.h Header.h
	g++ -c $(CFLAGS) $< -o $@

Normal.o: Normal.cpp Normal.h 
	g++ -c $(CFLAGS) $< -o $@

vMF.o: vMF.cpp vMF.h Header.h
	g++ -c $(CFLAGS) $< -o $@

FB4.o: FB4.cpp FB4.h Header.h
	g++ -c $(CFLAGS) $< -o $@

FB6.o: FB6.cpp FB6.h Header.h
	g++ -c $(CFLAGS) $< -o $@

Test.o: Test.cpp Test.h Header.h
	g++ -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *~ main gmon.out 

