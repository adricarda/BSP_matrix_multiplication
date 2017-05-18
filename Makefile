CC= gcc
CFLAGS= -Wall -O3 -g -std=c99 

LFLAGS= libmcbsp1.2.0.a -pthread -lrt -lm
#For MacOS use instead:
#LFLAGS= libmcbsp1.2.0.a -lm

OBJIP= bspinprod.o bspedupack.o
OBJBEN= bspbench.o bspedupack.o
OBJSORT= bspsort_test.o bspsort.o bspedupack.o
OBJODDEVEN = oddeven.o bspedupack.o
OBJMATRIXMUL = matrixMul.o bspedupack.o
OBJCANNON = cannon.o bspedupack.o

all: cannon matrixMul

inprod: $(OBJIP)
	$(CC) $(CFLAGS) -o inprod $(OBJIP) $(LFLAGS)

bench: $(OBJBEN)
	$(CC) $(CFLAGS) -o bench $(OBJBEN) $(LFLAGS)

sort: $(OBJSORT)
	$(CC) $(CFLAGS) -o sort $(OBJSORT) $(LFLAGS)

oddeven: $(OBJODDEVEN)
	$(CC) $(CFLAGS) -o oddeven $(OBJODDEVEN) $(LFLAGS)
	
matrixMul: $(OBJMATRIXMUL)
	$(CC) $(CFLAGS) -o matrixMul $(OBJMATRIXMUL) $(LFLAGS)

cannon: $(OBJCANNON)
	$(CC) $(CFLAGS) -o cannon $(OBJCANNON) $(LFLAGS)

bspinprod.o:  bspedupack.h
bspbench.o:   bspedupack.h
bspsort.o:    bspedupack.h
bspedupack.o: bspedupack.h
oddeven.o: bspedupack.h
matrixMul.o: bspedupack.h
cannon.o: bspedupack.h

.PHONY: clean
clean:
	rm -f $(OBJCANNON) $(OBJMATRIXMUL) $(OBJIP) $(OBJBEN) $(OBJSORT) $(OBJODDEVEN) inprod  bench sort oddeven matrixMul cannon
