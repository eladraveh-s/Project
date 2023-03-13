CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm
CC = gcc

spkmeans: spkmeans.o spkmeans.h
	$(CC) -o spkmeans spkmeans.o $(CFLAGS)

spkmeans.o: spkmeans.c
	$(CC) -c spkmeans.c $(CFLAGS)
