CC = gcc

mongrail: mongrail.c
	$(CC) -g -Wall -o mongrail mongrail.c -lm `pkg-config --cflags --libs glib-2.0`
clean:
	$(RM) mongrail.o
