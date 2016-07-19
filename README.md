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

First install the basics
```
apt-get install libglib2.0-dev cmake pkg-config gsl-bin libgsl-dev libnlopt-dev libsuitesparse-dev
```
Once those packages are installed its time to move on to `petsc`. We will configure `petsc` to install `MPICH` as our MPI installation. For the sake of example I will download everything to my home dir

```
cd ~
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.5.4.tar.gz .
tar zxvf petsc-3.5.4.tar.gz .
cd petsc-3.5.4/
./configure --download-fblaslapack --download-mpich
```
this will take a few minutes to run. next you are ready to make `petsc`
```
make PETSC_DIR=<yourhomedir>/petsc-3.5.4 PETSC_ARCH=arch-linux2-c-debug all
```
with `<yourhomedir>` replaced with the proper value of your home directory or where ever else you have downloaded petsc. once the compilation is complete, you can test it using `make test`. `petsc` relies on two environmental variables, `PETSC_DIR` and `PETSC_ARCH`. Its a good idea to add a couple lines to your .bash_profile such as
```
export PETSC_DIR=<yourhomedir>/petsc-3.5.4
export PETSC_ARCH=arch-linux2-c-debug
```
again inserting the proper value of `<yourhomedir>`

