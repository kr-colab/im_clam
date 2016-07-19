# im_clam

## Installation

`im_clam` is has several dependencies, thus installation will require install of multiple third party software packages.
Hopefully, with the help of a modern package manager like homebrew or apt-get this won't be too painful though. Here is a 
list of dependencies

- glib
- pkg-config
- csparse
- petsc (version 3.5)
- slepc (version 3.5)
- nlopt
- gsl

Please note that newer versions (>=3.6) of the petsc and slepc libraries break `im_clam` due to strange issues with backward compatibility. 