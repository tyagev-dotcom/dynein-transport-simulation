$ErrorActionPreference = "Stop"

New-Item -ItemType Directory -Force -Path "bin" | Out-Null

gcc -O3 -std=c11 -Wall -Wextra `
    "src/dynein_transport_simulation.c" `
    -lm -pthread `
    -o "bin/dynein_transport_simulation.exe"
