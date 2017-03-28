# peclet
A convection-diffusion solver written in C++ and based on deal.II

Author: Alexander G. Zimmerman <zimmerman@aices.rwth-aachen.de>

Doxygen generated HTML documentation: https://alexanderzimmerman.github.io/peclet/

# For developers:
## Versions

This is currently being tested with the following builds of deal.II:
- deal.II v8.4.2 built by candi (https://github.com/koecher/candi) on Ubuntu 14.04

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
