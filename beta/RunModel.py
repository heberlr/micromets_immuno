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
    callingModel = ['./melanoma', '/N/slate/hlimadar/melanoma_v3/output_S'+str("%06d"%Sample)+'_R'+str("%02d"%Replica)+'/config.xml']
    #callingModel = ['./melanoma', 'output/output_S'+str("%06d"%Sample)+'_R'+str("%02d"%Replica)+'/config.xml']
    # cache = subprocess.run( callingModel,universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cache = subprocess.run( callingModel,universal_newlines=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if ( cache.returncode != 0):
        print("Model output error! returned: "+ str(cache.returncode))
        os._exit(1)

if (len(sys.argv) != 3):
      print("Please provide 2 args: Initial and final sample (endpoint not included)")
      sys.exit(1)
initial_sample = int(sys.argv[1])
final_sample = int(sys.argv[2])

Samples = len(list(range(initial_sample,final_sample,1)))
numDataPerRank  = int(Samples/size)
mod = Samples%size
if ( mod != 0): numDataPerRank = numDataPerRank + 1
data = None
if rank == 0:
    data = np.array(list(range(initial_sample,final_sample,1)), dtype='d')
    if ( mod != 0):
      add = -1*np.ones((size-mod),dtype=int)
      data = np.concatenate((data,add),axis=None)
recvbuf = np.empty(numDataPerRank , dtype='d')
comm.Scatter(data, recvbuf, root=0)

for i in range(recvbuf.shape[0]):
    for idxRep in range(0,10,1):
        print('Rank: ',rank, ', recvbuf received: ',recvbuf[i], ', Replica: ', idxRep)
        if ( recvbuf[i] >= 0 ):
          model(recvbuf[i], idxRep)
