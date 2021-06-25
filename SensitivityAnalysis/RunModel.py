import os
import subprocess

from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def model(parameter=None, Sample=-1, Replica=-1): 

    cwd = os.getcwd()
    dir = "_calculation_" + parameter
    dir = os.path.join(cwd, dir)
    os.makedirs(dir)
    os.chdir(dir) 
    
    # Write input for simulation & execute
    callingModel = './ExecModel.exe '+str(rank)+' output/output_S'+str("%05d"%Sample)+'_R'+str("%03d"%Replica)+'.dat'
    print(callingModel)
    subprocess.run( callingModel)

    # After the program is finished, do something here with the output files
    # and return the data.
    data = str(parameter)+" S: " + str(Sample)+" R: " + str(Replica)
    
    os.chdir(cwd)
    
    
    return data

def scattering(samples, idSample, idReplica): #Scattering NumPy arrays
    sendbuf = None
    if rank == 0:
        sendbuf = np.empty(samples.shape, dtype='f')
        sendbuf = samples
    recvbuf = np.empty(samples.shape[1], dtype='f')
    comm.Scatter(sendbuf, recvbuf, root=0)   
    result = model(parameter=recvbuf,Sample=idSample, Replica=rank + (idReplica)*size)


if __name__ == '__main__':    
    # Reading samples
    NumReplicas = 12
    if (NumReplicas%size != 0): 
        print("Error: number of replicas needs to be divisible by number of processors! Replicas: " + str(NumReplicas)+" \# Processors: "+str(size))
    else:
        Samples = np.load('SA_samples.npy')
        for idSample in range(Samples.shape[0]):
            for idReplica in range(NumReplicas//size):
                TempSamples = np.array(size*[Samples[idSample]], dtype='f')
                scattering(TempSamples, idSample, idReplica)