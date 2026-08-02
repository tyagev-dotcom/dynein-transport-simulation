#!/usr/bin/env sh
set -eu

mkdir -p bin

cc="${CC:-gcc}"
cflags="${CFLAGS:--O3 -std=c11 -Wall -Wextra}"
ldflags="${LDFLAGS:--lm -pthread}"

"$cc" $cflags src/rigid_nanoparticle.c $ldflags -o bin/rigid_nanoparticle
"$cc" $cflags src/rigid_nanoparticle_low_persistence.c $ldflags -o bin/rigid_nanoparticle_low_persistence
"$cc" $cflags src/rigid_nanoparticle_flexible_input.c $ldflags -o bin/rigid_nanoparticle_flexible_input
"$cc" $cflags src/fluid_nanoparticle.c $ldflags -o bin/fluid_nanoparticle
"$cc" $cflags src/fluid_nanoparticle_optimized.c $ldflags -o bin/fluid_nanoparticle_optimized
