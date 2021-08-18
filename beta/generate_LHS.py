# This script provides parameter samples from latin hypercube. The script creates
# a new folder (subdirectory) for each set of parameters, makes changes to a default
# configuration (.xml) file using specified parameter values (in an accompanying .txt file),
# copies the new config file into the new folder

import xml.etree.ElementTree as ET
from shutil import copyfile
from scipy.stats import qmc
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def generate_parSamples(parameters, default_value, variation, Samples_number,Replicas_number,constrain,fileOut,seed=12345):
    # Resampling samples2 (simple random sampling) to satifisfy the constrain samples1 < samples2
    def ResamplingWithConstrain(samples1, samples2, b1, h1, b2, h2):
        indx = np.argwhere(samples1 >= samples2).flatten()
        if (len(indx) > 0):
            samples2[indx] = samples1[indx] + (np.random.uniform(0,1,len(indx))*(h2-samples1[indx]))
        if ( (samples1 < samples2).all() ):
            print(str(len(indx)) +" resampling done using constrained simple random sampling!")

    file = open(fileOut, "w")
    # Latin Hypercube Sampling with multi-dimensional uniformity ( (N,M) N variables with M samples each)
    sampler = qmc.LatinHypercube(d=len(parameters),seed=seed)
    samples = sampler.random(n=Samples_number)
    min = default_value - variation*default_value
    max = default_value + variation*default_value

    #Transform unit hypercube in real values
    for j in range(0, len(parameters)):
        samples[:,j] = min[j] + samples[:,j]*(max[j]-min[j])

    #Resampling with constrain x1 < x2
    for j in range(0, len(constrain)):
        idx1 = constrain[j][0]
        idx2 = constrain[j][1]
        ResamplingWithConstrain(samples[:,idx1],samples[:,idx2],min[idx1],max[idx1],min[idx2],max[idx2])

    #Write file with samples
    for sample_id in range(Samples_number):
        for replica_id in range(Replicas_number):
            folder = 'output_S'+str("%06d"%sample_id)+'_R'+str("%02d"%replica_id)
            file.write("folder"+" "+folder+"\n")
            # Set of parameters
            for id_par in range(0, len(parameters)):
                file.write(parameters[id_par]+" "+str(samples[sample_id,id_par])+"\n")
            file.write("custom_save_time true\n")
            file.write("bitmap_save true\n")
            file.write("#"+"\n")
    file.close()

def generate_configXML(params_file):
    xml_file_in = 'config/PhysiCell_settings.xml'
    xml_file_out = 'config/tmp.xml'
    copyfile(xml_file_in,xml_file_out)
    tree = ET.parse(xml_file_out)
    xml_root = tree.getroot()
    first_time = True
    output_dirs = []
    with open(params_file) as f:
        for line in f:
            print(len(line),line)
            if (line[0] == '#'):
                continue
            (key, val) = line.split()
            if (key == 'folder'):
                if first_time:  # we've read the 1st 'folder'
                    first_time = False
                else:  # we've read  additional 'folder's
                    # write the config file to the previous folder (output) dir and start a simulation
                    print('---write (previous) config file and start its sim')
                    tree.write(xml_file_out)
                xml_file_out = val + '/config.xml'  # copy config file into the output dir
                output_dirs.append(val)
            if ('.' in key):
                k = key.split('.')
                uep = xml_root
                for idx in range(len(k)):
                    uep = uep.find('.//' + k[idx])  # unique entry point (uep) into xml
                uep.text = val
            else:
                if (key == 'folder' and not os.path.exists(val)):
                    print('creating ' + val)
                    os.makedirs(val)
                xml_root.find('.//' + key).text = val

    tree.write(xml_file_out)
    print(output_dirs)

def plot_samples(filename, label1, label2, default_value1, default_value2):
    file = open(filename)
    n_line = 0
    value1 = []
    value2 = []
    for line in file:
        n_line += 1
        line = line.strip().split(' ')
        if label1 in line:
            #print(str(n_line)+" "+line[0]+" "+line[1])
            value1.append(float(line[1]))
        if label2 in line:
            #print(str(n_line)+" "+line[0]+" "+line[1])
            value2.append(float(line[1]))
    file.close()
    plt.scatter(np.asarray(value1),np.asarray(value2),s=1)
    plt.axhline(y=default_value2, color='r', linestyle='--',alpha=0.5)
    plt.axvline(x=default_value1, color='r', linestyle='--',alpha=0.5)
    plt.xlabel(label1)
    plt.ylabel(label2)
    plt.show()

if __name__ == '__main__':
    parameters = np.array(["macrophage_max_recruitment_rate","macrophage_recruitment_min_signal","macrophage_recruitment_saturation_signal","DC_max_recruitment_rate","DC_recruitment_min_signal","DC_recruitment_saturation_signal","DC_leave_prob","TC_death_rate","T_Cell_Recruitment","DM_decay"])
    default_value = np.array([4e-9, 0.1, 0.3, 2e-9, 0.1, 0.3, 3.3e-9, 1.4e-6, 1.1e-4, 3.5e-4])
    variation = 1.0
    constrain = ([1,2],[4,5]) # parameters[1] < parameters[2] and parameters[4] < parameters[5]
    file = "ParameterSamples1.txt"
    Samples_number = 10000
    Replicas_number = 10
    # Generate samples from Latin Hypercube
    generate_parSamples(parameters, default_value, variation, Samples_number,Replicas_number,constrain, file)
    # Create .xml and folder to each simulation
    ##generate_configXML(file)
    # Plot samples
    plot_samples(file, parameters[1], parameters[2], default_value[1], default_value[2])
    plot_samples(file, parameters[4], parameters[5], default_value[4], default_value[5])
    plot_samples(file, parameters[0], parameters[3], default_value[0], default_value[3])
    plot_samples(file, parameters[6], parameters[7], default_value[6], default_value[7])
    plot_samples(file, parameters[8], parameters[9], default_value[8], default_value[9])
