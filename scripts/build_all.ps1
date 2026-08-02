$ErrorActionPreference = "Stop"

New-Item -ItemType Directory -Force -Path "bin" | Out-Null

$targets = @(
    @{ Source = "src/rigid_nanoparticle.c"; Output = "bin/rigid_nanoparticle.exe" },
    @{ Source = "src/rigid_nanoparticle_low_persistence.c"; Output = "bin/rigid_nanoparticle_low_persistence.exe" },
    @{ Source = "src/rigid_nanoparticle_flexible_input.c"; Output = "bin/rigid_nanoparticle_flexible_input.exe" },
    @{ Source = "src/fluid_nanoparticle.c"; Output = "bin/fluid_nanoparticle.exe" },
    @{ Source = "src/fluid_nanoparticle_optimized.c"; Output = "bin/fluid_nanoparticle_optimized.exe" }
)

foreach ($target in $targets) {
    Write-Host "Building $($target.Output)"
    gcc -O3 -std=c11 -Wall -Wextra $target.Source -lm -pthread -o $target.Output
}
