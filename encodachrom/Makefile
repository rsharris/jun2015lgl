CFLAGS = -O3 -Wall -Wextra -Werror
CC     = gcc
LDLIBS = -lm

PROGS= encodachrom decodachrom

all: ${PROGS}

encodachrom:   encodachrom.o
decodachrom:   decodachrom.o

clean:
	rm -f ${PROGS}
	rm -f core *.o

