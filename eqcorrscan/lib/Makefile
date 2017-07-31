# Makefile for normalise utils C routines - not used in build, only for testing

CC = gcc
CFLAGS = -pthread -Wsign-compare -DNDEBUG -g -fwrapv -O3 -Wall -fPIC
FFTWFLAGS = -lfftw3 -lfftw3_threads
FFTWLOC = `pkg-config --cflags --libs-only-L fftw3`

all: multi_corr

multi_corr: multi_corr.c
	$(CC) $(CFLAGS) $(FFTWLOC) $(FFTWFLAGS) -fopenmp -g -O0 multi_corr.c -o multi_corr -lm -lgomp

clean:
	rm -f norm
	rm -f multi_corr