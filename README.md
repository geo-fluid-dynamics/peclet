# peclet
an unsteady scalar convection-diffusion solver based on the finite element C++ library deal.II

[![Build Status](https://travis-ci.org/alexanderzimmerman/peclet.svg?branch=master)](https://travis-ci.org/alexanderzimmerman/peclet) (<b>Continuous integration status</b>; click the button to go to Travis-CI)

Author: Alexander G. Zimmerman <zimmerman@aices.rwth-aachen.de>

Doxygen generated HTML documentation: https://alexanderzimmerman.github.io/peclet/

# For users:
## Run pre-built version on docker image

Pull the image from https://hub.docker.com/r/zimmerman/peclet/

    docker pull zimmerman/peclet

# For developers:
## Versions

This is currently being tested with the following builds of deal.II:
- deal.II v8.5.pre from docker image dealii/dealii:v8.5.pre.4-gcc-mpi-fulldepsmanual-debugrelease (as shown in peclet/Dockerfile)

## Design notes
The Peclet class is implemented entirely with header files. This reduces the structural complexity of the code and can increase programming productivity, but it leads to longer compile times. A header-only approach would be impractical for the deal.II library itself; but in this small project's experience, the header-only approach is more than adequate. Most notably, this simplifies working with C++ templates.

## Build

    git clone git@github.com:alexanderzimmerman/peclet.git

    mkdir build

    cd build

    cmake ../peclet

    make test
    
## Documentation
The Doxygen generated HTML docs are hosted in the standard GitHub fashion on the gh-pages branch.

For initial set up with your local clone, follow the procedure found in https://github.com/m-a-d-n-e-s-s/madness/issues/104

Then whenever commiting to the master branch, update the gh-pages branch as follows:

    git push origin master

    doxygen

    cd doc/html

    git add *

    git commit -m "Refreshed HTML doc"

    git push origin gh-pages
