# loga

Loga is an experimental tool for extracting log message templates and placeholders from raw log files, 
without requiring structured or labeled training data or extraction rules. It performs 
a pipeline of sequence-alignment–based operations to align log messages and the patterns derived from them.

> Initial code extracted and refactored from the [`neel/prova`](https://github.com/neel/prova) repository.


## Usage

loga is a command line application. Only necessary input is the log file. 
It is assumed that the logfile is a sequence of ASCII messages seperated by new line characters.
A terminal that supports ANSI colors is necessary to view the output.

```bash
./loga Logfile.log
```

![Apache Log Demo](https://github.com/user-attachments/assets/baec328f-f0e8-41d1-9701-8ad3c59f44bb)


## Building

This project can be built in two ways:

1. **Using system packages only** (recommended on Linux if your distro provides everything).
2. **Using [vcpkg](https://github.com/microsoft/vcpkg)** to fetch all C/C++ dependencies.

The CMake build itself is the same in both cases; the only difference is how you install the libraries.

---

### Dependencies

Common requirements: CMake, C++20 compiler 

C/C++ libraries:

- Boost (`icl`, `graph`, `asio`, `program_options`)
- Armadillo
- cereal
- igraph
---

### Build with system packages (no vcpkg)

Install dependencies:

##### Arch Linux

```bash
sudo pacman -Syu
sudo pacman -S base-devel cmake git boost cereal igraph
yay -S armadillo
```

#### Compile 

```bash
git clone https://github.com/neel/loga.git
cd loga
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

### Build using vcpkg

```bash
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.sh   # Linux or Mac
.\bootstrap-vcpkg.bat  # Windows
```

#### Compile 

```bash
git clone https://github.com/neel/loga.git
cd loga

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../../vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build .
```

#### Compile (Windows without MinGW) 

```bash
cmake -S .. -B . -DCMAKE_BUILD_TYPE=Release "-DCMAKE_TOOLCHAIN_FILE=../../vcpkg/scripts/buildsystems/vcpkg.cmake" -DVCPKG_TARGET_TRIPLET=x64-windows-static-md -DVCPKG_HOST_TRIPLET=x64-windows
cmake --build build --config Release --parallel
```