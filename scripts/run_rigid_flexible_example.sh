#!/usr/bin/env sh
set -eu

cp examples/Input-NP-flex.txt Input-NP-flex.txt

if [ ! -x bin/rigid_nanoparticle_flexible_input ]; then
    ./scripts/build_all.sh
fi

./bin/rigid_nanoparticle_flexible_input
