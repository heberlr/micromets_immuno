import os
import subprocess

from mpi4py import MPI
import numpy as np
import sys

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def model(Sample=-1, Replica=-1):
    # Write input for simulation & execute
    callingModel = ['./micromets_lung', 'PatientVariation/output_S'+str("%06d"%Sample)+'_R'+str("%02d"%Replica)+'/config.xml']
    # callingModel = ['./micromets_lung', '/N/slate/hlimadar/micromets_lung_v3/output_S'+str("%06d"%Sample)+'_R'+str("%02d"%Replica)+'/config.xml']
    # callingModel = ['./micromets_lung', 'output/output_S'+str("%06d"%Sample)+'_R'+str("%02d"%Replica)+'/config.xml']
    # cache = subprocess.run( callingModel,universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cache = subprocess.run( callingModel,universal_newlines=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if ( cache.returncode != 0):
        print("Model output error! returned: "+ str(cache.returncode))
        os._exit(1)

if __name__ == '__main__':
    if (len(sys.argv) < 4):
        print("Error: Number of args!")
        print("option 1) [Initial ID sample] [Final ID sample] sequential")
        print("option 2) [ ID_sample_1  ID_sample_2 ... ID_sample_n]")
        sys.exit(1)

    if ( sys.argv[3] == 'sequential' ):
        initial_sample = int(sys.argv[1])
        final_sample = int(sys.argv[2])
        Samples = np.array(list(range(initial_sample,final_sample+1,1)), dtype='d')
    else:
        Samples = np.zeros(shape=(len(sys.argv)-1), dtype='d')
        for sample_index in range(len(Samples)):
            Samples[sample_index] = int(sys.argv[sample_index+1])

    NumReplicates = 10
    numDataPerRank  = int(len(Samples)/size)
    mod = len(Samples)%size
    if ( mod != 0): numDataPerRank = numDataPerRank + 1
    data = None

    if rank == 0:
        data = np.array(Samples, dtype='d')
        if ( mod != 0):
          add = -1*np.ones((size-mod),dtype='d')
          data = np.concatenate((data,add),axis=None)

    recvbuf = np.empty(numDataPerRank , dtype='d')
    comm.Scatter(data, recvbuf, root=0)
    for i in range(recvbuf.shape[0]):
        for idxRep in range(0,NumReplicates,1):
            print('Rank: ',rank, ', recvbuf received: ',recvbuf[i], ', Replica: ', idxRep)
            if ( recvbuf[i] >= 0 ):
              model(recvbuf[i], idxRep)
