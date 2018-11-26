
CC = gcc
CFLAGS = -O3 -Wall -Wno-unused-result

LIBS = -lm -lgsl -lgslcblas

HEADERS = funcs.h
OBJECTS = funcs.o


%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: Task1

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

