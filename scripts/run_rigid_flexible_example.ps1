$ErrorActionPreference = "Stop"

Copy-Item -Force "examples/Input-NP-flex.txt" "Input-NP-flex.txt"

if (!(Test-Path "bin/rigid_nanoparticle_flexible_input.exe")) {
    & ".\scripts\build_all.ps1"
}

& ".\bin\rigid_nanoparticle_flexible_input.exe"
