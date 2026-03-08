# Introduction-to-Image-Understanding-2
This repository contains the solution for Homework 2 of CENG 391 - Introduction to Image Understanding.

## Assignment Description
This project extends the `ceng391::Image` library built in the previous assignment by introducing spatial filtering, image derivatives, geometric transformations, and PNG file support.

## Key Features Implemented

### 1. Image Filtering (Smoothing & Box Filters)
* **Gaussian Smoothing:** Implemented 1D spatial Gaussian smoothing functions (`smooth_x`, `smooth_y`) and a combined 2D smoothing function (`smooth`) using dynamically generated Gaussian kernels.
* **Box Filtering:** Implemented separable box filters (`box_filter_x`, `box_filter_y`, and `box_filter`) for efficient image blurring.

### 2. Image Derivatives
* Implemented functions to compute the spatial derivatives of the image in the X and Y directions (`deriv_x`, `deriv_y`) using standard 1D derivative kernels (e.g., `[-1, 0, 1]` and `[1, 2, 1]`). 

### 3. Image Rotation & Interpolation
* Implemented an image `rotate` function that allows rotating the image by a specified angle (`theta`).
* **Interpolation Methods:** Supports both **Nearest Neighbor** and **Bilinear** interpolation to handle sub-pixel mappings during rotation.
* **Center Rotation:** Includes a flag to rotate the image around its center rather than the origin.

### 4. PNG Support integration
* Transitioned from basic PNM format to **PNG format** using `libpng`.
* Implemented `xload_png` and `xsave_png` for robust reading and writing of grayscale, RGB, and RGBA PNG images.

## Prerequisites
* **CMake** (v3.10+)
* **C++17** compatible compiler
* **libpng** (Required for PNG loading/saving)

## How to Build and Run
1. Clone this repository to your local machine.
2. Open your terminal and navigate to the project directory.
3. Create a build directory, compile, and build the project using CMake:
   ```bash
   mkdir build
   cd build
   cmake ..
   make
