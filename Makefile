CC ?= gcc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra
LDFLAGS ?= -lm -pthread

BIN_DIR := bin
TARGET := $(BIN_DIR)/dynein_transport_simulation
SOURCE := src/dynein_transport_simulation.c

.PHONY: all clean

all: $(TARGET)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(TARGET): $(SOURCE) | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@

clean:
	rm -rf $(BIN_DIR)
