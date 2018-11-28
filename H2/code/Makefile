
CC = gcc
CFLAGS = -O3 -Wall

LIBS = -lm -lgsl -lgslcblas

HEADERS = funcs.h
OBJECTS = funcs.o


%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: Task1 Task2

Task1: $(OBJECTS) main_T1.c
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Task2: $(OBJECTS) main_T2.c
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

# $(PROGRAMS): $(OBJECTS) main_T1.c
#	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

