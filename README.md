# im_clam

## Installation

`im_clam` is has several dependencies, thus installation will require install of multiple third party software packages.
Hopefully, with the help of a modern package manager like homebrew or apt-get this won't be too painful though. Here is a 
list of dependencies

- CMake
- atlas (or another BLAS and LAPACK package)
- glib
- pkg-config
- SuiteSparse
- petsc (version 3.5)
- slepc (version 3.5)
- nlopt
- gsl
- a version of MPI 

Please note that newer versions (>=3.6) of the petsc and slepc libraries break `im_clam` due to strange issues with backward compatibility. One day I will get around to changing all this code, but it is not a current priority. One thing that eases the process of installing dependencies is that `petsc` will take care of many of the packages if you tell it to when you configure the install. I'll demonstrate this below.

## Install on Linux
I focus on linux here as most compute clusters suitable for running `im_clam` will, I assume, be using some flavor of linux. Here I go through the steps necessary on a clean Ubuntu 16.04 install. 
```
sudo apt-get install libglib2.0-dev cmake 
```