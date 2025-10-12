# image-processing

* A small project that combines my implementations of tasks from the Computer Graphics course, taught in the 3rd year at
  university

## Requirements

* **OS:** Linux
* **Build Tools:**
    * A C++23 compliant compiler (GCC 13+ or Clang 16+)
    * `CMake` (version 3.23+)
* **Dependencies:**
    * All dependencies are managed via `vcpkg` and are installed automatically during the build process

## Building the Project

1. **Clone the repository with submodules:**
   ```bash
   git clone --recursive https://github.com/antilopinae/image-processing
   cd image-processing
   ```

2. **Configure and build the project:**
   The `build.sh` script automates this process
   ```bash
   ./build.sh
   cd build
   ```
   Alternatively, you can run the commands manually:
   ```bash
   # Configure the project, specifying the static triplet
   cmake -B build -S .

   # Build the project
   cmake --build build -j$(nproc)
   cd build
   ```

## Usage

### General format:

```bash
./image-processing <LAB_COMMAND> <SUBCOMMAND> [OPTIONS...]
```

---

### Lab 1 — `lab1`

#### 1. Building a round grayscale image

Creates a PNG image of the specified size with a gray circle in the middle

```bash
./image-processing lab1 circle-gray -o <output.png>
```

**Arguments:**

* `-o`, `--output` — output filepath to PNG-file

*By default, the image is generated in 512×512 size with a circle radius of 0.45 of the width*

Example of usage:

```bash
./image-processing lab1 circle-gray -o circle.png
```

---

#### 2. Blending two images using the alpha channel

Performs byte-by-byte mixing of two 8 bpp PNG images, using the third image as the alpha channel

```bash
./image-processing lab1 blend \
  --input-first <A.png> \
  --input-second <B.png> \
  --input-alpha <alpha.png> \
  -o <output.png>
```

**Arguments:**

* `--input-first` — first input PNG-image
* `--input-second` — second input PNG-image
* `--input-alpha` — input image with alpha (PNG, 8 bpp)
* `-o`, `--output` — output filepath to PNG-file

---

### Lab 2 — `lab2`

(unimplemented)

```bash
./image-processing lab2 --input <input.png> -o <output.png>
```

---

### Help

Print help:

```bash
./image-processing --help
```

Print lab1:

```bash
./image-processing lab1 --help
```

Print lab1 blend:

```bash
./image-processing lab1 blend --help
```

And others...