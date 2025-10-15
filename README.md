# image-processing

* A small project that combines my implementations of tasks from the Computer Graphics course, taught in the 3rd year at
  university

## Requirements

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
  --image-first <A.png> \
  --image-second <B.png> \
  --image-alpha <alpha.png> \
  -o <output.png>
```

**Arguments:**

* `--image-first` — first input PNG-image
* `--image-second` — second input PNG-image
* `--image-alpha` — input image with alpha (PNG, 8 bpp)
* `-o`, `--output` — output filepath to PNG-file

---

### Lab 2 — `lab2`

Converting an 8 bpp image to n bpp (n < 8) using the Floyd-Stenberg error scattering algorithm

```bash
./image-processing lab2 --input <input.png> -o <output.png> -n <number_of_levels>
```

**Arguments:**

* `--input` — path to input PNG-image for scattering
* `-n`, `--n-levels` — amount of bpp to convert
* `-o`, `--output` — output filepath to PNG-file

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

---

# Results

### Laboratory 1: Print gray circle

![example_lab1_circle](./assets/ex-circle.png)

### Laboratory 1: Blending two images

#### First image

![example_lab1_blend_image1](./assets/ex-blend-image1.png)

#### Second image

![example_lab1_blend_image2](./assets/ex-blend-image2.png)

#### Alpha channel

![example_lab1_blend_alpha](./assets/ex-blend-alpha.png)

#### Result

![example_lab1_blend_result](./assets/ex-blend-result.png)

---

### Laboratory 2

#### Input image

![example_lab2_floyd](./assets/ex-floyd.png)

#### Output with 2 bpp

![example_lab2_floyd_result](./assets/ex-floyd2.png)

#### Output with 3 bpp

![example_lab2_floyd_result](./assets/ex-floyd3.png)

#### Output with 4 bpp

![example_lab2_floyd_result](./assets/ex-floyd4.png)

#### Output with 8 bpp

![example_lab2_floyd_result](./assets/ex-floyd8.png)
