gfortran -O3 -g -fbacktrace -cpp -Wall -Wno-maybe-uninitialized -Wno-unused-dummy-argument -fcheck=all \
../src/lists.f90  \
../src/strings.f90 \
../src/stop.f90 \
../src/parameters.f90 \
../src/arrayutilities.f90 \
../src/boundaries.f90 \
../src/scalarboundaries.f90 \
interpolate_new_grid.f90 \
-o ~/bin/interpolate_new_grid \
"$@"
