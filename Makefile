CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -lm
CC = gcc

spkmeans: spkmeans.o spkmeansmodule.o spkmeans.h
	$(CC) -o spkmeans spkmeans.o spkmeansmodule.o $(CFLAGS)

spkmeans.o: spkmeans.c
	$(CC) -c spkmeans.c $(CFLAGS)

spkmeansmodule.o: spkmeansmodule.c
	$(CC) -c spkmeansmodule.c $(CFLAGS)
