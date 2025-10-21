CC=gcc
CFLAGS=-O2 -Wall -Wextra
LDFLAGS=-lm

SRC=main.c tsp_io.c tsp_distance.c tsp_matrice.c tsp_force_brute.c tsp_force_brute_matrice.c
OBJ=$(SRC:.c=.o)
BIN=tsp

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(BIN)

