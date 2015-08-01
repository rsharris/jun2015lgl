CFLAGS = -O3 -Wall -Wextra -Werror
CC     = gcc
LDLIBS = -lm

default:   chrom_avg

chrom_avg: chrom_avg.o utilities.o

incFiles = utilities.h

%.o: %.c Makefile ${incFiles}
	${CC} -c ${CFLAGS} $< -o $@

clean: cleano
	rm -f chrom_avg

cleano:
	rm -f *.o
