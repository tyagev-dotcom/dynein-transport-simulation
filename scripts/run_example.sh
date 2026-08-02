#!/usr/bin/env sh
set -eu

if [ ! -x bin/dynein_transport_simulation ]; then
    ./scripts/build_all.sh
fi

./bin/dynein_transport_simulation
