# Makefile for normalise utils C routines - not used in build, only for testing

CC = g++
CFLAGS = -pthread -Wsign-compare -DNDEBUG -g -fwrapv -O3 -Wall -fPIC
FFTWFLAGS = -lfftw3f
FFTWLOC = `pkg-config --cflags --libs-only-L fftw3`

all: norm multi_corr

norm: norm.c
	$(CC) $(CFLAGS) -g norm.c -o norm -lm

multi_corr: multi_corr.cpp
	$(CC) $(CFLAGS) $(FFTWLOC) $(FFTWFLAGS) -fopenmp -g multi_corr.cpp -o multi_corr -lm

clean:
	rm -f norm
	rm -f multi_corr