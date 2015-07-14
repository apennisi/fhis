# FHIS
Fast Hsv Image Segmentation (FHIS) Library is an OpenCV based C++ adaptation of the original 
Matlab code designed for performing an accurate segmentation in real-time. FHIS creates a simple 
representation of the image by using the Delaunay triangulation. A HSV threshold is used in order 
to find similar triangles, obtaining an accurate segmented image. FHIS exploits OpenCV 2 and CGAL functions. 
The library is based on the work [1] realized by Camillo Taylor and Anthony Cowley of University of Pennsylvania.

## Requirements

FHIS requires the following packeges to build:

* OpenCV (< 3.0)
* CGal

## How to build

AT works under Linux and Mac Os environments. We recommend a so-called out of source build 
which can be achieved by the following command sequence:

* mkdir build
* cd build
* cmake ../
* make -j\<number-of-cores+1\>

## How to use

Once the build phase has been successfully, you can use FGHISby launching the _fhis_
file, located in the _bin_ folder.

For more information, you can visit the following link: [here](http://www.dis.uniroma1.it/~pennisi/FHIS-Image_Segmentation_Library.html).
