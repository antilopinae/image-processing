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

### Lab3 - `lab3`

Building and filling polygons.

Implement the following functions:

1. Drawing straight line segments with a thickness of 1 pixel.
2. Displaying the polygon on the screen.
3. Definitions of the polygon type: simple or complex (i.e. with self-intersections), convex or non-convex.
4. Filling the polygon using the even-odd and non-zero-winding rules for determining
   whether a pixel belongs to a polygon.

### Lab4 - `lab4`

Draw Bezier curves, clipping straight line segments with a convex polygon.

Implement the following functions:

1. Construction of Bezier curves of the third order.
2. Clipping straight line segments with a convex polygon using the algorithm
   Kirusa-Beka.

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

### Laboratory 3

#### Output images

![example_lab3_polygon_edges](./assets/polygon_edges.png)

![example_lab3_polygon_even_odd](./assets/polygon_filled1.png)

![example_lab3_polygon_non_zero](./assets/polygon_filled2.png)

![example_lab3_star_even_odd](./assets/lab3_star_EO.png)

![example_lab3_star_non_zero](./assets/lab3_star_NZW.png)

![example_lab3_lines](./assets/lab3_thick_lines.png)

### Laboratory 4

#### Output images

#### Cross-segments

![example_lab4_cross_square](./assets/segment_cross_square.png)

![example_lab4_inside_square](./assets/segment_inside_square.png)

![example_lab4_outside_square](./assets/segment_outside_square.png)

![example_lab4_cyrus_beck](./assets/lab4_cyrus_beck.png)

![example_lab4_sutherland_hondgman](./assets/lab4_sutherland_hodgman.png)

#### Bezier curves

![example_lab4_bezier_basic](./assets/bezier_basic.png)

![example_lab4_bezier_linearity](./assets/bezier_linearity.png)

![example_lab4_bezier_symmetric](./assets/bezier_symmetric.png)

![example_lab4_bezier_adaptive](./assets/lab4_bezier_adaptive.png)

### Laboratory 5

#### Output images

##### Parallel

![example_lab5_parallepiped](./assets/animation/lab5_parallel_000.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_001.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_002.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_003.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_004.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_005.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_006.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_007.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_008.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_009.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_010.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_011.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_012.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_013.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_014.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_015.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_016.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_017.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_018.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_019.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_020.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_021.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_022.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_023.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_024.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_025.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_026.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_027.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_028.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_029.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_030.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_031.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_032.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_033.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_034.png)
![example_lab5_parallepiped](./assets/animation/lab5_parallel_035.png)

##### Perspective

![example_lab5_parallepiped](./assets/animation/lab5_perspective_000.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_001.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_002.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_003.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_004.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_005.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_006.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_007.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_008.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_009.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_010.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_011.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_012.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_013.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_014.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_015.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_016.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_017.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_018.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_019.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_020.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_021.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_022.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_023.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_024.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_025.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_026.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_027.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_028.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_029.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_030.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_031.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_032.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_033.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_034.png)
![example_lab5_parallepiped](./assets/animation/lab5_perspective_035.png)

### Homework 1

#### Output images

Original image

![example_hw1_original](./assets/hw1_original.png)

Result

![example_hw1](./assets/hw1.png)

![example_hw1](./assets/hw1_complex_star.png)

### Homework 2

#### Output images

![example_hw2](./assets/hw2.png)

![example_hw2](./assets/hw2_arc_gallery.png)

### Homework 3

#### Input images

![example_hw3_input](./assets/ex-floyd.png)

#### Output images

![example_hw3](./assets/hw3.png)