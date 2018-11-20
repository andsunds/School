
CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm -lgsl -lgslcblas

HEADERS = initfcc.h alpotential.h funcs.h
OBJECTS = initfcc.o alpotential.o funcs.o MD_main.o
PROGRAM = MD

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

