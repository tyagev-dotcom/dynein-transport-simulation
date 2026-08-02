$ErrorActionPreference = "Stop"

if (!(Test-Path "bin/dynein_transport_simulation.exe")) {
    & ".\scripts\build_all.ps1"
}

& ".\bin\dynein_transport_simulation.exe"
