CC      = gcc
LD      = gcc
CFLAGS  = -ggdb3 -Wall -fPIC -I./include
LDFLAGS =

all: matrix.so ops.so

clean:
	rm *.so *.o

# make shared library out of the object file
%.so: %.o
	$(LD) -shared $(LDFLAGS) -o $@ $< utils.o
	rm utils.o

# compile source file to object file
%.o: src/%.c
	$(CC) $(CFLAGS) -fPIC -c src/utils.c $< 
