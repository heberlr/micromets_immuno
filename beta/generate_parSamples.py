# This script provides parameter samples from latin hypercube
import sys
import lhsmdu
import numpy as np

if (len(sys.argv) < 3):
  usage_str = "Usage: %s <Samples_number> <Replicas_number>" % (sys.argv[0])
  print(usage_str)
  print("e.g.:  python generate_parSamples.py 10000 10")
  exit(1)
else:
   Samples_number = int(sys.argv[1])
   Replicas_number = int(sys.argv[2])

file = open("beta/ParameterSamples.txt", "w")



parameters = ["macrophage_max_recruitment_rate","macrophage_recruitment_min_signal","macrophage_recruitment_saturation_signal","DC_max_recruitment_rate","DC_recruitment_min_signal","DC_recruitment_saturation_signal","DC_leave_prob","TC_death_rate","T_Cell_Recruitment","DM_decay"]
default_value = [4e-9, 0.1, 0.3, 2e-9, 0.1, 0.3, 3.3e-9, 1.4e-6, 1.1e-4, 3.5e-4]

# Latin Hypercube Sampling with multi-dimensional uniformity ( (N,M) N variables with M samples each)
samples = lhsmdu.sample(len(parameters), Samples_number )
variation = 1.0
min = np.asarray(default_value) - variation*np.asarray(default_value)
max = np.asarray(default_value) + variation*np.asarray(default_value)

#Transform unit hypercube in real values
for j in range(0, len(parameters)):
    samples[j,:] = min[j] + samples[j,:]*(max[j]-min[j])

#Constrain:  macrophage_recruitment_min_signal < macrophage_recruitment_saturation_signal
if (np.greater(samples[1,:], samples[2,:]).any()):
    indexes = np.argwhere(np.greater(samples[1,:], samples[2,:])).flatten()
    temp = samples[1,indexes]
    samples[1,indexes] = samples[2,indexes]
    samples[2,indexes] = temp
#Constrain: DC_recruitment_min_signal < DC_recruitment_saturation_signal
if (np.greater(samples[4,:], samples[5,:]).any()):
    indexes = np.argwhere(np.greater(samples[4,:], samples[5,:])).flatten()
    temp = samples[4,indexes]
    samples[4,indexes] = samples[5,indexes]
    samples[5,indexes] = temp

# Write file with samples
for sample_id in range(Samples_number):
    for replica_id in range(Replicas_number):
        folder = 'output_S'+str("%06d"%sample_id)+'_R'+str("%02d"%replica_id)
        file.write("folder"+" "+folder+"\n")
        # Set of parameters
        for id_par in range(0, len(parameters)):
            file.write(parameters[id_par]+" "+str(samples[id_par,sample_id])+"\n")
        file.write("custom_save_time true\n")
        file.write("#"+"\n")
file.close()
