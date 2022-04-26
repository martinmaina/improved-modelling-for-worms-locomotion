# Locomotion in C.elegans 

Code for running a simluation of a C. elegans integrated neuro-mechanical locomotion model written in `C++`.

If using the mechanical model please cite: Cohen, Netta and Ranner, Thomas (2017), A new computational method for a model of C. elegans biomechanics: Insights into elasticity and locomotion performance.  arXiv e-prints 1702.04988, https://arxiv.org/abs/1702.04988

If using the neuromechanical model please cite: Denham, Jack E. and Ranner, Thomas and Cohen, Netta (2018), Intrinsic and extrinsic modulation of C. elegans locomotion,  bioarxiv e-prints 10.1101/312132,    https://www.biorxiv.org/content/early/2018/05/01/312132

## Prerequisites

`C++` compiler with support for `C++11` (e.g. GNU C++ Compiler version 4.8.5 or above), `CMake` version 2.8 or above,
GNU `make`,
Basic Linear Algebra Subroutines (BLAS) (e.g. atlas-blas),
Suitesparse

(See more detailed description below).

## Synopsis

  ./worm-model ../data/model-parameters [options]

     Run simluation with parameters from `model-parameters`
     file. Extra options can be used which overwrite model-parameters
     using `--parameter=new-value`. Outputs are:

     * worm_*.vtu all solution variables output in a format readable by Paraview
     * Activation.txt - neural activation
     * Beta.txt  - bending moment applied to body
     * Kappa.txt  - scalar body curvature
     * feedback.txt - feedback measure by neurons
     * pointx - cartesian x-coordinate of all points along body
     * pointy - cartesian y-coordinate of all points along body
     * energy.txt

   The above output files contain N+1 columns and dt(end - transient) rows.
   Parameters N (#meshpoints), dt (integration time step), end (simulation end time) 
   and transient (cut off initial transient) are specified the the model-parameters file.

## Authors
Admin

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.


## Extra details for building prerequisites

### Ubuntu Steps

Install dependencies:
```
apt-get update
apt-get install -y g++
apt-get install -y libatlas-base-dev
apt-get install -y libsuitesparse-dev
apt-get install -y git
```

Clone and build repository:
```
git clone https://github.com/martinmaina/improved-modelling-for-worms-locomotion
cd improved-modelling-for-worms-locomotion
mkdir build && cd build
cmake -D CMAKE_BUILD_TYPE=Release ..
make
```

Optionally test:
```
ctest
```

### Mac steps (untested)

Install dependencies

1. git clone <<repository URL>>
2. brew install cmake
3. chmod -R 777 python folder (TODO: Hack - this shouldn't be necessary)
4. brew install gcc
5. remove rt from CMakeLists.txt
6. sudo port -v install SuiteSparse
7. sudo port install opencv

## Running

1. paraview &
2. Open output/worm...vtu
3. Play

## Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

## Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
