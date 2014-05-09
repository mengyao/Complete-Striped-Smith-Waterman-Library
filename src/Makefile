CC = gcc
CFLAGS := -Wall -Werror -O3 -pipe #-pg
LOBJS = ssw.o
PROG = ssw_test
EXAMPLE = example

.PHONY: all clean

all: $(PROG) $(EXAMPLE)

$(PROG): main.c

$(EXAMPLE): example.c

$(PROG) $(EXAMPLE): $(LOBJS)
	$(CC) -o $@ $^ $(CFLAGS) -lm -lz

ssw.o: ssw.c ssw.h
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	-rm -f $(LOBJS) $(PROG) $(EXAMPLE) *~
