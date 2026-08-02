#!/usr/bin/env sh
set -eu

mkdir -p bin

cc="${CC:-gcc}"
cflags="${CFLAGS:--O3 -std=c11 -Wall -Wextra}"
ldflags="${LDFLAGS:--lm -pthread}"

"$cc" $cflags src/dynein_transport_simulation.c $ldflags -o bin/dynein_transport_simulation
