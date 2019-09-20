# peclet
an unsteady scalar convection-diffusion solver based on the finite element C++ library deal.II

[![Build Status](https://travis-ci.org/alexanderzimmerman/peclet.svg?branch=master)](https://travis-ci.org/alexanderzimmerman/peclet) (<b>Continuous integration status</b>; click the button to go to Travis-CI)

Author: Alexander G. Zimmerman <zimmerman@aices.rwth-aachen.de>

Doxygen generated HTML documentation: https://alexanderzimmerman.github.io/peclet/

This code supports [the author's masters thesis](https://arxiv.org/abs/1909.08882).

# For users:
## Run pre-built version on docker image
Get the free community edition of Docker here: https://www.docker.com/community-edition

Pull the image from https://hub.docker.com/r/zimmerman/peclet/ and run the container with docker

    docker run -ti zimmerman/peclet:latest
    
Or run the container with access to a shared folder (shared between the host and the container)

    docker run -ti -v $(pwd):/home/dealii/shared zimmerman/peclet:latest
    
If you plan to use this container repeatedly, then instead use this command to also give it a name

    docker run -ti -v $(pwd):/home/dealii/shared --name peclet zimmerman/peclet:latest

After exiting the container, you can start it again with

    docker start peclet
    
You can confirm that the container is running with

    docker ps
    
or list all containers (running or not) with

    docker ps -a

To enter a bash terminal inside of the running container

    docker start peclet
    
    docker exec -ti -u dealii peclet /bin/bash -l

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

The procedure for keeping the HTML docs updated is rough. To initially create the gh-pages branch, we followed the outline in an [issue](https://github.com/m-a-d-n-e-s-s/madness/issues/104) from another repository, which uses the ideas from [here](http://rickfoosusa.blogspot.de/2011/10/howto-use-doxygen-with-github.html) and [here](https://gist.github.com/chrisjacob/825950). Since any of these links may break, in short the procedure was

    cd peclet
    
    mkdir doc
    
    mkdir doc/html
    
    cd doc/html
    
    git clone git@github.com:alexanderzimmerman/peclet.git .
    
    git checkout -b gh-pages
    
    git branch -d master
    
    git rm -r *

    git commit "Removed everything from gh-pages branch"
    
    cd ../..
    
    git checkout master
    
    doxygen
    
    cd html/doc
    
    git checkout gh-pages
    
    git add *
    
    git commit "Added all HTML documentation"
    
    git push origin gh-pages

You can skip many of those steps for initial set up with your local clone. Simply make the target documents directory and clone the gh-pages branch inside of it.

    cd peclet

    mkdir doc
    
    mkdir doc/html
    
    cd doc/html
    
    git clone git@github.com:alexanderzimmerman/peclet.git .
    
    git checkout gh-pages

Then whenever commiting to the master branch, update the gh-pages branch as follows:

    cd peclet
    
    git push origin master

    doxygen

    cd doc/html

    git add *

    git commit -m "Refreshed HTML doc"

    git push origin gh-pages
