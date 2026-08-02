CC ?= gcc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra
LDFLAGS ?= -lm -pthread

BIN_DIR := bin

TARGETS := \
	$(BIN_DIR)/rigid_nanoparticle \
	$(BIN_DIR)/rigid_nanoparticle_low_persistence \
	$(BIN_DIR)/rigid_nanoparticle_flexible_input \
	$(BIN_DIR)/fluid_nanoparticle \
	$(BIN_DIR)/fluid_nanoparticle_optimized

.PHONY: all clean rigid rigid-lp rigid-flex fluid fluid-opt

all: $(TARGETS)

rigid: $(BIN_DIR)/rigid_nanoparticle
rigid-lp: $(BIN_DIR)/rigid_nanoparticle_low_persistence
rigid-flex: $(BIN_DIR)/rigid_nanoparticle_flexible_input
fluid: $(BIN_DIR)/fluid_nanoparticle
fluid-opt: $(BIN_DIR)/fluid_nanoparticle_optimized

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/rigid_nanoparticle: src/rigid_nanoparticle.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@

$(BIN_DIR)/rigid_nanoparticle_low_persistence: src/rigid_nanoparticle_low_persistence.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@

$(BIN_DIR)/rigid_nanoparticle_flexible_input: src/rigid_nanoparticle_flexible_input.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@

$(BIN_DIR)/fluid_nanoparticle: src/fluid_nanoparticle.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@

$(BIN_DIR)/fluid_nanoparticle_optimized: src/fluid_nanoparticle_optimized.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(LDFLAGS) -o $@

clean:
	rm -rf $(BIN_DIR)
