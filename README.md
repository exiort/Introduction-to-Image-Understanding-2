# Introduction-to-Image-Understanding-2

This repository contains the solution for Homework 1 of CENG 391 - Introduction to Image Understanding.

# Assignment Description
The assignment consists of three main exercises. The purpose of each exercise is to extend the functionality of the ceng391::Image class. Here's a brief overview of what each exercise requires:

# Exercise 1 - Color Support for the ceng391::Image Class
Added support for specifying the number of color channels per pixel (1 for grayscale, 3 for RGB, or 4 for RGBA) in the constructor.
Implemented RGB and RGBA versions for the set and set_rect functions to work with multi-channel data.
Modified the write_pnm function to save RGB data in the binary PPM format.
# Exercise 2 - Loading PNM Images
Added a new member function Image::read_pnm to read image contents from PGM or PPM binary formats.
When reading color images, it creates a four-channel RGBA image with alpha values initialized to 255.
# Exercise 3 - Color Conversion
Implemented three new member functions: Image::to_rgb, Image::to_grayscale, and Image::to_rgba.
These functions convert the image contents from the current number of channels to the desired format.
Includes formulas for converting between RGB/RGBA and grayscale, taking care of underflow and overflow.
# How to Use the Code
* Clone this repository to your local machine.
* Open the project in your preferred development environment.
* Compile and build the project.
* Use the provided functions in the ceng391::Image class to manipulate and convert images.
