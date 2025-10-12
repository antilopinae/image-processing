# image-processing

* A small project that combines my implementations of tasks from the Computer Graphics course, taught in the 3rd year at
  university.

## Requirements

* **OS:** Linux
* **Build Tools:**
    * A C++23 compliant compiler (GCC 13+ or Clang 16+).
    * `CMake` (version 3.23+).
* **Dependencies:**
    * All dependencies are managed via `vcpkg` and are installed automatically during the build process.

## Building the Project

1. **Clone the repository with submodules:**
   ```bash
   git clone --recursive <your_repository_url>
   cd <project_folder_name>
   ```

2. **Configure and build the project:**
   The `build.sh` script automates this process.
   ```bash
   ./build.sh
   ```
   Alternatively, you can run the commands manually:
   ```bash
   # Configure the project, specifying the static triplet
   cmake -B build -S .

   # Build the project
   cmake --build build -j$(nproc)
   ```

## Usage

```bash
sudo ./build/image-processing
```
