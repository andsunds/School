CC = gcc # c- compiler
CFLAGS = -O3 # compiler flags such as -O3 -W
LIBS = -lm # libraries such as -lm -lgsl -lgslcblas
OBJECTS = harm_osc.o main.o # object files to produce
PROGRAM = a # name of executable to be built

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

# compile program
all: $(PROGRAM)
$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)


clean:
ifneq ("$(wildcard *.o)","")
	rm $(OBJECTS)
endif
ifneq ("$(wildcard $(PROGRAM))","")
	rm $(PROGRAM)
endif
