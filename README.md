# loga

Code seperated from repo neel/prova

## Usage

loga is a command line application. Only necessary input is the log file. 
It is assumed that the logfile is a sequence of ASCII messages seperated by new line characters.
A terminal that supports ANSI colors is necessary to view the output.

```bash
./loga Logfile.log
```



## Building

This project can be built in two ways:

1. **Using system packages only** (recommended on Linux if your distro provides everything).
2. **Using [vcpkg](https://github.com/microsoft/vcpkg)** to fetch all C/C++ dependencies.

The CMake build itself is the same in both cases; the only difference is how you install the libraries.

---

### Dependencies

Common requirements:

- CMake
- A C++20 compiler 
- Git 

C/C++ libraries:

- Boost (`thread`, `program_options`)
- Armadillo
- cereal
- igraph
---

### Build with system packages (no vcpkg)

Install dependencies:

##### Ubuntu 

```bash
sudo apt update
sudo apt install build-essential cmake git libboost-thread-dev libboost-program-options-dev libarmadillo-dev libcereal-dev libigraph-dev
```

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

./bootstrap-vcpkg.sh   # Linux or Mac
.\bootstrap-vcpkg.bat  # Windows
```

#### Compile 

```bash
git clone https://github.com/neel/loga.git
cd loga

mkdir build
cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE=../../vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build .
```
