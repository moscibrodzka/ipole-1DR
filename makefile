CC = gcc
LDFLAGS = -lm -lgsl -lgslcblas 

SRC = \
polar_main.c polar_lib.c

OBJ = \
polar_main.o polar_lib.o

all: $(OBJ) makefile
	 $(CC) $(CFLAGS) -o ipole-1DR $(OBJ) $(LDFLAGS) 

$(OBJ) : makefile decs.h

clean:
	rm *.o
cleanup:
	rm ipole-1DR



