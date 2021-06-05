

# Compile c++ program
g++ .\c-OpenMO.cpp -fopenmp -o ExecModel

# Runnign python script
numprocs=4
mpiexec -n numprocs python -m mpi4py RunExternal.py