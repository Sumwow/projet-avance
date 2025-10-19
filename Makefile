CC      = gcc
CFLAGS  = -O2 -Wall -Wextra
LDFLAGS = -lm

SRC = main.c tsp_io.c tsp_distance.c tsp_matrice.c
OBJ = $(SRC:.c=.o)

tsp: $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJ) tsp

