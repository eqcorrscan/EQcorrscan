# Makefile for normalise utils C routines - not used in build, only for testing

CC = gcc
CFLAGS = -pthread -Wsign-compare -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -lm

all: norm

norm: norm.c
	$(CC) $(CFLAGS) -g norm.c -o norm