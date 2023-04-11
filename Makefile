CC = cc
VPATH = src include
FLAGS = 

all: gendiplo mongrail

debug: FLAGS += -g -Wall
debug: all

gendiplo: gendiplo.c
	$(CC) $(FLAGS) -o gendiplo $< -I include -lm

mongrail: mongrail.c
	$(CC) $(FLAGS) -o mongrail $< -lm `pkg-config --cflags --libs glib-2.0`


clean:
	$(RM) mongrail gendiplo
