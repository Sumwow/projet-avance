CC = gcc
CFLAGS = -O2 -Wall -Wextra
LDFLAGS = -lm

# Tous les fichiers sources C
SRC := $(wildcard *.c)

# Sources communes (tous sauf les deux main)
SRC_COMMON := $(filter-out main.c main_dbg.c, $(SRC))

# Exécutables
BIN = tsp
BIN_DBG = tsp_dbg

# Objets correspondants
OBJ := $(SRC_COMMON:.c=.o) main.o
OBJ_DBG := $(SRC_COMMON:.c=.o) main_dbg.o

all: $(BIN) $(BIN_DBG)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

$(BIN_DBG): $(OBJ_DBG)
	$(CC) $(OBJ_DBG) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(OBJ_DBG) $(BIN) $(BIN_DBG)

