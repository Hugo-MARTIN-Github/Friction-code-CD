all: matrices.so

matrices.so: matrices.c
		$(CC) -Wall -g -std=c99 -fPIC -shared -o $@ $? -lc

clean:
		rm -f matrices.so
