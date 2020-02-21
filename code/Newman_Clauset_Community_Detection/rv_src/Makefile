CFLAGS = -O
CC = gcc
LIBS = -lm -lgsl

metadata:bp-dc-groups.o readwg2.o
	$(CC) $(CFLAGS) -o metadata $(LIBS) bp-dc-groups.o readwg2.o

metadata.o:bp-dc-groups.c readwg2.h network.h Makefile
	$(CC) $(CFLAGS) -c bp-dc-groups.c

readwg2.o:readwg2.c network.h Makefile
	$(CC) $(CFLAGS) -c readwg2.c
