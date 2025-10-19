.PHONY: all clean

all: tsp

tsp:
	@gcc -O2 -Wall -Wextra -o tsp main.c tsp_io.c tsp_distance.c -lm

clean:
	@rm -f tsp
