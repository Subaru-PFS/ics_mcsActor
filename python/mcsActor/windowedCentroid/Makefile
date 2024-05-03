TARGET = test_centroid
LIBS = -lm -lcfitsio -L.
CC = gcc
CFLAGS = -Wall -fPIC -g

.PHONY: clean all default

default: $(TARGET)
all: default

SOURCES = centroid_win.c newCentroid.c test_centroid.c fitswrappers.c
OBJECTS = $(SOURCES:.c=.o)
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o
