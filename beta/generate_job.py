import sys
import os
import subprocess
import numpy as np
import xml.etree.ElementTree as ET
from shutil import copyfile
import pathlib

def create_xml_file(xml_file_in, xml_file_out, parameters):
    if not os.path.exists(xml_file_out):
        print(pathlib.Path(xml_file_out).parent)
        os.makedirs(pathlib.Path(xml_file_out).parent)
    copyfile(xml_file_in,xml_file_out)
    tree = ET.parse(xml_file_out)
    xml_root = tree.getroot()
    for key in parameters.keys():
        val = parameters[key]
        print(key,val)
        if ('.' in key):
            k = key.split('.')
            uep = xml_root
            for idx in range(len(k)):
                uep = uep.find('.//' + k[idx])  # unique entry point (uep) into xml
                # print(k[idx])
            uep.text = val
        else:
            if (key == 'folder' and not os.path.exists(val)):
                print('creating ' + val)
                os.makedirs(val)

            xml_root.find('.//' + key).text = val
    tree.write(xml_file_out)

def createFolderConfig(Samples):
    macrophage_max_recruitment_rate = [0,4e-9,8e-9]
    DC_max_recruitment_rate = [0,2e-9,4e-9]
    NumReplicates = 10
    for idSample in Samples:
        for idDC,DC_value in enumerate(DC_max_recruitment_rate):
            for idM,M_value in enumerate(macrophage_max_recruitment_rate):
                for idReplicate in range(NumReplicates):
                    folderName = "PatientVariation/output_S"+str("%06d"%(idSample+1e5*(idDC*len(DC_max_recruitment_rate) + idM + 1)))+'_R'+str("%02d"%idReplicate)
                    parameters = {'folder': folderName,'DC_max_recruitment_rate':str(DC_value), 'macrophage_max_recruitment_rate':str(M_value), 'only_population': 'true'}
                    # create_xml_file('../CPDT/PatientsNC/output_S'+str("%06d"%idSample)+'_R00/config.xml', folderName+'/config.xml', parameters)
                    create_xml_file('PatientsAnalysis/output_S'+str("%06d"%idSample)+'_R00/config.xml', folderName+'/config.xml', parameters)

def create_fileSamples(ID_Job,Samples,NumNodes,NumTaskPerNode,NumCPUsPerTask):
  original_stdout = sys.stdout
  sys.stdout = open('job_'+str("%06d"%ID_Job)+'.sh','w')

  print ("#!/bin/bash\n")
  print ("#SBATCH --mail-user=hlimadar@iu.edu")
  print ("#SBATCH --job-name=CPDTs"+str("%06d"%ID_Job))
  print ("#SBATCH -p general")
  print ("#SBATCH -o CPDTs_%j.txt")
  print ("#SBATCH -e CPDTs_%j.err")
  print ("#SBATCH --nodes="+str(NumNodes))
  print ("#SBATCH --ntasks-per-node="+str(NumTaskPerNode))
  print ("#SBATCH --cpus-per-task="+str(NumCPUsPerTask))
  print ("#SBATCH --time=1-0:00:00")
  print ("#SBATCH --mail-type=FAIL")
  print ("#SBATCH --mem=16G\n")
  print ("module load python/3.6.11")
  print ("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")

  run = 'python /beta/RunModel.py'
  for idSample in Samples:
    run += ' '+str(int(idSample))
  print ("srun --cpu-bind=sockets "+run)

  sys.stdout = original_stdout

def create_fileIndiv(ID_Job,idSample, idRep, NumCPUsPerTask):
  original_stdout = sys.stdout
  sys.stdout = open('job_'+str("%06d"%ID_Job)+'.sh','w')

  print ("#!/bin/bash\n")
  print ("#SBATCH --mail-user=hlimadar@iu.edu")
  print ("#SBATCH --job-name=CPDTs"+str("%06d"%ID_Job))
  print ("#SBATCH -p general")
  print ("#SBATCH -o CPDTs_%j.txt")
  print ("#SBATCH -e CPDTs_%j.err")
  print ("#SBATCH --nodes=1")
  print ("#SBATCH --cpus-per-task="+str(NumCPUsPerTask))
  print ("#SBATCH --time=0-1:00:00")
  #print ("#SBATCH --mail-type=FAIL,BEGIN,END")
  print ("#SBATCH --mem=16G\n")
  print ("module load python/3.6.11")
  print ("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")

  run = './PulmonaryMicroMets /N/slate/hlimadar/PulmonaryMicroMets_v3/output_S'+str("%06d"%idSample)+'_R'+str("%02d"%idRep)+'/config.xml'
  print ("srun --cpu-bind=sockets "+run)

  sys.stdout = original_stdout


def sub_job(jobname):
  command = ['sbatch',jobname]
  cache = subprocess.Popen( command)
  #os.remove(jobname)

def CreatePatientsVariation():
    createFolderConfig([515, 723, 1810, 2539, 2723, 2958, 3160, 5379, 6040, 7642]) # (No control)
    createFolderConfig([93, 1420, 2602, 3081, 4884, 4972, 5633, 6868, 7956, 8222]) # (Marginal control)


if __name__ == '__main__':
  CreatePatientsVariation()
  exit(0)
  NumNodes = 2
  NumTaskPerNode = 3
  NumCPUsPerTask = 8
  # Samples
  SamplesByRank = 17
  NumRanks = NumNodes*NumTaskPerNode

  Samples = np.array([100515., 100723., 101810., 102539., 102723., 102958., 103160.,
       105379., 106040., 107642., 100093., 101420., 102602., 103081.,
       104884., 104972., 105633., 106868., 107956., 108222., 200515.,
       200723., 201810., 202539., 202723., 202958., 203160., 205379.,
       206040., 207642., 200093., 201420., 202602., 203081., 204884.,
       204972., 205633., 206868., 207956., 208222., 300515., 300723.,
       301810., 302539., 302723., 302958., 303160., 305379., 306040.,
       307642., 300093., 301420., 302602., 303081., 304884., 304972.,
       305633., 306868., 307956., 308222., 400515., 400723., 401810.,
       402539., 402723., 402958., 403160., 405379., 406040., 407642.,
       400093., 401420., 402602., 403081., 404884., 404972., 405633.,
       406868., 407956., 408222., 500515., 500723., 501810., 502539.,
       502723., 502958., 503160., 505379., 506040., 507642., 500093.,
       501420., 502602., 503081., 504884., 504972., 505633., 506868.,
       507956., 508222., 600515., 600723., 601810., 602539., 602723.,
       602958., 603160., 605379., 606040., 607642., 600093., 601420.,
       602602., 603081., 604884., 604972., 605633., 606868., 607956.,
       608222., 700515., 700723., 701810., 702539., 702723., 702958.,
       703160., 705379., 706040., 707642., 700093., 701420., 702602.,
       703081., 704884., 704972., 705633., 706868., 707956., 708222.,
       800515., 800723., 801810., 802539., 802723., 802958., 803160.,
       805379., 806040., 807642., 800093., 801420., 802602., 803081.,
       804884., 804972., 805633., 806868., 807956., 808222., 900515.,
       900723., 901810., 902539., 902723., 902958., 903160., 905379.,
       906040., 907642., 900093., 901420., 902602., 903081., 904884.,
       904972., 905633., 906868., 907956., 908222.])
  NumJobs = np.ceil(len(Samples)/(NumRanks*SamplesByRank))
  SamplesJobs = np.array_split(Samples, NumJobs)
  for i in range(len(SamplesJobs)):
    create_fileSamples(i,SamplesJobs[i],NumNodes,NumTaskPerNode,NumCPUsPerTask )
  exit(0)

  if (len(sys.argv) < 4):
    print("Error: Number of args!")
    print("option 1) [ID sample] [ID replicate] individual")
    print("option 2) [Initial ID sample] [Final ID sample] sequential")
    print("option 3) [ ID_sample_1  ID_sample_2 ... ID_sample_n]")
    sys.exit(1)
  if ( len(sys.argv) == 4 ):
    if ( sys.argv[3] == 'individual' ):
    	idSample = int(sys.argv[1])
    	idRep= int(sys.argv[2])
    	ID_job = idSample*10 + idRep
    	create_fileIndiv(ID_job,idSample,idRep,NumCPUsPerTask)
    	sub_job('job_'+str("%06d"%ID_job)+'.sh')
    else:
        if (  sys.argv[3] == 'sequential' ):
            Samples = np.array(list(range(int(sys.argv[1]),int(sys.argv[2])+1,1)), dtype='d')
        else:
            Samples = np.zeros(shape=(len(sys.argv)-1), dtype='int')
            for sample_index in range(len(Samples)):
                Samples[sample_index] = int(sys.argv[sample_index+1])
        NumJobs = np.ceil(len(Samples)/(NumRanks*SamplesByRank))
        SamplesJobs = np.array_split(Samples, NumJobs)
        for i in range(len(SamplesJobs)):
          create_fileSamples(i,SamplesJobs[i],NumNodes,NumTaskPerNode,NumCPUsPerTask )
          # sub_job('job_'+str("%06d"%i)+'.sh') # Submit job
