import sys
import os
import subprocess

def create_fileSamples(ID_Job,idSampleInitial,idSampleFinal):
  original_stdout = sys.stdout 
  sys.stdout = open('job_'+str("%06d"%ID_Job)+'.sh','w')
  
  print ("#!/bin/bash\n")
  print ("#SBATCH --mail-user=hlimadar@iu.edu")
  print ("#SBATCH --job-name=CPDTs"+str("%06d"%ID_Job))
  print ("#SBATCH -p general")
  print ("#SBATCH -o CPDTs_%j.txt")
  print ("#SBATCH -e CPDTs_%j.err")
  print ("#SBATCH --nodes=10")
  print ("#SBATCH --ntasks-per-node=3")
  print ("#SBATCH --cpus-per-task=8")
  print ("#SBATCH --time=1-0:00:00")
  print ("#SBATCH --mail-type=FAIL")
  print ("#SBATCH --mem=16G\n")
  print ("module load python/3.6.11")
  print ("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
  
  run = 'python RunModel.py '+str(idSampleInitial)+' '+str(idSampleFinal)
  print ("srun --cpu-bind=sockets "+run)

  sys.stdout = original_stdout 

def create_fileIndiv(ID_Job,idSample, idRep):
  original_stdout = sys.stdout 
  sys.stdout = open('job_'+str("%06d"%ID_Job)+'.sh','w')
  
  print ("#!/bin/bash\n")
  print ("#SBATCH --mail-user=hlimadar@iu.edu")
  print ("#SBATCH --job-name=CPDTs"+str("%06d"%ID_Job))
  print ("#SBATCH -p general")
  print ("#SBATCH -o CPDTs_%j.txt")
  print ("#SBATCH -e CPDTs_%j.err")
  print ("#SBATCH --nodes=1")
  print ("#SBATCH --cpus-per-task=8")
  print ("#SBATCH --time=0-1:00:00")
  #print ("#SBATCH --mail-type=FAIL,BEGIN,END")
  print ("#SBATCH --mem=16G\n")
  print ("module load python/3.6.11")
  print ("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
  
  run = './melanoma /N/slate/hlimadar/melanoma_v3/output_S'+str("%06d"%idSample)+'_R'+str("%02d"%idRep)+'/config.xml'
  print ("srun --cpu-bind=sockets "+run)

  sys.stdout = original_stdout 


def sub_job(jobname):
  command = ['sbatch',jobname]
  cache = subprocess.Popen( command)
  #os.remove(jobname)

if __name__ == '__main__':
  if ( len(sys.argv) == 4 ):
    if ( sys.argv[3] == 'individual' ):
    	idSample = int(sys.argv[1])
    	idRep= int(sys.argv[2])
    	ID_job = idSample*10 + idRep
    	create_fileIndiv(ID_job,idSample,idRep)
    	sub_job('job_'+str("%06d"%ID_job)+'.sh')
    if (  sys.argv[3] == 'samples' ):
        NumJobs= int(sys.argv[1]) # 20
        SamplesByRank = int(sys.argv[2])# 17
        NumRanks = 30
        for i in range(NumJobs):
          InitialSample = i*SamplesByRank*NumRanks
          FinalSample = InitialSample + SamplesByRank*NumRanks
          if ( FinalSample > 10000): FinalSample = 10000
          create_fileSamples(i,InitialSample ,FinalSample )
          sub_job('job_'+str("%06d"%i)+'.sh') # Submit job    
  else:
    print("Error number of args!")
    print("option 1) [ID sample] [ID replicate] individual")
    print("option 2) [Number of jobs] [Samples by rank] samples")
    exit(0)
