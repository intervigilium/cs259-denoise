# Makefile

TARGET=riciandenoise3d
SRC= riciandenoise3d.c
PAPIDIR=/mnt/jc5/CS259/papi

CC=gcc

INCFLAGS= -I$(PAPIDIR)
LDFLAGS= -L/mnt/jc5/CS259/papi/ -lm -lutil_papi -lpapi
CFLAGS= -pg -g

default: $(TARGET)

clean:
	rm -f $(TARGET)

$(TARGET):
	$(CC) $(CFLAGS) $(INCFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)
