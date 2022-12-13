mongrail: read_indv_data.c
	gcc -g -Wall -o mongrail read_indv_data.c -lm `pkg-config --cflags --libs glib-2.0`
