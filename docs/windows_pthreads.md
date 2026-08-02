# Windows pthreads Setup Notes

The simulation source uses POSIX threads through `pthread.h`.

## Option 1: MSYS2 / MinGW-w64

1. Install MSYS2.
2. Install a MinGW-w64 GCC toolchain with pthread support.
3. Open the MSYS2 MinGW shell.
4. From the repository root, run:

```sh
make
```

## Option 2: Visual Studio with pthreads-win32

1. Create a C console project in Visual Studio.
2. Add one desired source file from `src/`.
3. Add pthreads-win32 headers to the C/C++ include path.
4. Add the pthreads-win32 library directory to the linker library path.
5. Add the pthreads library, such as `pthreadVC2.lib`, to linker input dependencies.
6. Place the matching runtime DLL, such as `pthreadVC2.dll`, beside the executable.

Keep third-party pthreads binaries outside the repository unless their license and redistribution terms have been reviewed.
