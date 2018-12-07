
CC = gcc
CFLAGS = -O3 -Wall 

LIBS = -lm -lgsl -lgslcblas

HEADERS = funcs.h 
OBJECTS = funcs.o main.o

PROGRAM = run

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM) 

$(PROGRAM): $(OBJECTS) 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

