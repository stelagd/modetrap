# modetrap
This application shows the effects of mode trapping on the eigenfrequency spacings of a vibrating string.

## Requirements

This package requires a fortran compiler, *f2py* to convert fortran into a loadable python module, matplotlib, pylab, and the TkAgg backend. It has been developed and tested on OSX 10.2.2 using the [Anaconda python distribution](https://store.continuum.io/cshop/anaconda/), with the following configuration:
- matplotlib 1.4.3
- numpy 1.9.2
- python 2.7.9
- gfortran: GNU Fortran (GCC) 4.9.3 20150202
- f2py 2.0

To build the modetrap binary file, type

    make modetrap_sub

This creates the file *modetrap_sub.so*. 

## Running

To start the application type 

    python modetrap_slider.py 1

for the one bead case, or 

    python modetrap_slider.py 2

for the two bead case. Click on the cases and adjust the sliders to reproduce the frequency spacing patterns shown, commonly referred to as _mode trapping_.
