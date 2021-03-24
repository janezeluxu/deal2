export TMPDIR=/tmp
DYLD_LIBRARY_PATH=
CC=/Users/zeluxu/mpi/bin/mpicc
DEAL_II_WITH_MPI=On
#cmake -DCMAKE_C_COMPILER="/Users/zeluxu/mpi/bin/mpicc" -DDEAL_II_DIR=../../../dealii-master/build/ .
cmake -DCMAKE_C_COMPILER="/Users/zeluxu/mpi/bin/mpicc" -DCMAKE_BUILD_TYPE="Release" -DDEAL_II_DIR=/Applications/deal.II.app/Contents/Resources/ .
make 
