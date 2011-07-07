CC=			gcc
CFLAGS=		-g -Wall -O3 #-m64 #-arch ppc
LOBJS=		ssw.o	
PROG=		ssw_test
all:$(PROG)

.PHONY:all clean cleanlocal
ssw_test:$(LOBJS) main.c 
		$(CC) $(CFLAGS) main.c -o $@ $(LOBJS) -lm -lz
ssw.o:ssw.h
cleanlocal:
		rm -fr *.o $(PROG) *~ 

clean:cleanlocal-recur


