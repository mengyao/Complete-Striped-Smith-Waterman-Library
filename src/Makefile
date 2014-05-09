CC = gcc
CXX = g++
CFLAGS := -Wall -Werror -O3 -pipe #-pg
CXXFLAGS := $(CFLAGS)
LOBJS = ssw.o
LCPPOBJS = ssw_cpp.o
PROG = ssw_test
EXAMPLE = example_c
EXAMPLE_CPP = example_cpp

.PHONY: all clean

all: $(PROG) $(EXAMPLE) $(EXAMPLE_CPP)

$(PROG): main.c kseq.h

$(EXAMPLE): example.c

$(PROG) $(EXAMPLE): $(LOBJS)
	$(CC) -o $@ $(filter-out %.h,$^) $(CFLAGS) -lm -lz

$(EXAMPLE_CPP): example.cpp $(LOBJS) $(LCPPOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) -lm -lz

ssw.o: ssw.c ssw.h
	$(CC) -c -o $@ $< $(CFLAGS)

ssw_cpp.o: ssw_cpp.cpp ssw_cpp.h ssw.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	-rm -f $(LOBJS) $(LCPPOBJS) $(PROG) $(EXAMPLE) $(EXAMPLE_CPP) *~
