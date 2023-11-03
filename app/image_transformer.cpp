#include <cstdlib>
#include <iostream>

#include "ceng391/image.hpp"

using namespace std;
using ceng391::Image;

int main(int argc, char** argv)
{
    if (argc != 2) {   
        cerr << "Usage: " << argv[0] << " <filename.png> " << endl;
        exit(EXIT_FAILURE);
    }

    Image *img = Image::from_png(argv[1], false);
    Image *rotated = Image::from_png(argv[1], false);

    img->smooth_x(1.5f);
    img->xsave_png("/tmp/smooth_x.png");
    
    img = Image::from_png(argv[1], false);
    img->smooth_y(1.5f);
    img->xsave_png("/tmp/smooth_y.png");

    img = Image::from_png(argv[1], false);
    img->smooth(1.5f, 1.5f);
    img->xsave_png("/tmp/smooth.png");

    img = Image::from_png(argv[1], false);
    img->rotate(0.25f, rotated, false, false);
    rotated->xsave_png("/tmp/nearest_neighbor_1.png");

    img = Image::from_png(argv[1], false);
    img->rotate(0.25f, rotated, false, true);
    rotated->xsave_png("/tmp/nearest_neighbor_2.png");

    img = Image::from_png(argv[1], false);
    img->rotate(0.25f, rotated, true, false);
    rotated->xsave_png("/tmp/bilinear_1.png");

    img = Image::from_png(argv[1], false);
    img->rotate(0.25f, rotated, true, true);
    rotated->xsave_png("/tmp/bilinear_2.png");

    delete img;
    delete rotated;
    
    return EXIT_SUCCESS;
}