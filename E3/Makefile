
CC = gcc
CFLAGS = -O3 -Wall -Wno-unused-result

LIBS = -lm -lgsl -lgslcblas

HEADERS = initfcc.h alpotential.h funcs.h fft_func.h
OBJECTS = initfcc.o alpotential.o funcs.o fft_func.o


%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: Task1 Task2 Task3 main_Prod.c

Task1: $(OBJECTS) main_T1.c
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Task2: $(OBJECTS) main_T2.c
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Task3: $(OBJECTS) main_T3.c
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Prod: $(OBJECTS) main_Prod.c
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

# $(PROGRAMS): $(OBJECTS) main_T1.c
#	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

