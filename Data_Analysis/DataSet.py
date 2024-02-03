import numpy as np
import pandas as pd
import pickle
import wget, os
import torchvision.transforms as transforms
from scipy.stats import ttest_ind, mannwhitneyu
from fitter import Fitter

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from matplotlib_venn import venn3_unweighted
import seaborn as sns

NameParameters = ["macrophage_max_recruitment_rate","macrophage_recruitment_min_signal","macrophage_recruitment_saturation_signal","DC_max_recruitment_rate","DC_recruitment_min_signal","DC_recruitment_saturation_signal","DC_leave_prob","TC_death_rate","T_Cell_Recruitment","DM_decay"]
NC_Patients_samples = np.array([515, 723, 1810, 2539, 2723, 2958, 3160, 5379, 6040, 7642]) # (No control)
MC_Patients_samples = np.array([93, 1420, 2602, 3081, 4884, 4972, 5633, 6868, 7956, 8222]) # (Marginal control)

def loading_pkl_file(filename):
    try:
        with open(filename, "rb") as f:
            return pickle.load(f)
    except Exception as ex:
        print("Error during unpickling object (Possibly unsupported):", ex)


def download_file(local_path, link):
    path, filename = os.path.split(local_path)
    if not os.path.exists(local_path):
        print('Downloading from %s, this may take a while...' % link)
        wget.download(link)
        print()

def Loading_dataset():
    # Download dataset frames
    local_path = 'DataFrames.pickle'
    link = 'https://figshare.com/ndownloader/files/38135118'
    download_file(local_path, link)
    df = loading_pkl_file(local_path)
    return df['input'], df['output']

def Loading_AditionalDataset():
    # Download the analysis of 20 patients
    local_path = 'dfPatientVariation.pickle'
    link = 'https://figshare.com/ndownloader/files/44381864'
    download_file(local_path, link)
    df = loading_pkl_file(local_path)
    return df

def Loading_UMAP_traj():
    # Download the UMAP of trajectories
    local_path = 'Clust_UMAP_15_01_2.pickle' # trajectories of all simulations
    link = 'https://figshare.com/ndownloader/files/38481611'
    download_file(local_path, link)
    # return loading_pkl_file(local_path)
    return pd.read_pickle(local_path)

def Loading_UMAP_avg():
    local_path2 = 'Patient_UMAP_15_01_2.pickle' # trajectory averages of all simulations
    link2 = 'https://figshare.com/ndownloader/files/38898579'
    download_file(local_path2, link2)
    return loading_pkl_file(local_path2)
    
def ZscoreImagesToTensor(dataframe): # Normalize images (zscore) and convert to tensor (pytorch pattern)
    images = np.asarray(dataframe['imagens'].tolist(), dtype="float64") # shape (1e4, 5, 5, 3)
    images_norm = np.zeros(shape = (images.shape[0],images.shape[3],images.shape[1],images.shape[2]), dtype="float64") # shape (1e4, 3, 5, 5)
    mean = np.mean(images,axis=(0,1,2)) # mean channels
    std = np.std(images,axis=(0,1,2)) # std channels
    # define custom transform function
    transform_norm = transforms.Compose([transforms.ToTensor(), transforms.Normalize(mean, std)])
    for idSample in range(images.shape[0]): images_norm[idSample,:,:,:] = transform_norm(images[idSample,:,:,:])
    return images_norm

def Classifier_Simulations(dataframe):
    Label = []
    thrs_LowNumCancerCell = 100 # threshold to define marginal and significant control
    thrs_HighNumCancerCell = 1500 # threshold to define marginal and no control
    df = dataframe.loc[dataframe['simulation_time'] ==  14400.0] # last frame
    images_norm = ZscoreImagesToTensor(df)
    listImagensNorm = [images_norm[i,:,:,:] for i in range(images_norm.shape[0])]
    LiveCancerCells = df['live_cancer'].to_numpy()
    for idx in range(len(df)):
        if ( LiveCancerCells[idx] >=  thrs_HighNumCancerCell ): Label.append('NC') # no control
        if ( (LiveCancerCells[idx] >=  thrs_LowNumCancerCell) & (LiveCancerCells[idx] <  thrs_HighNumCancerCell) ): Label.append('MC') # marginal control
        if ( LiveCancerCells[idx] <  thrs_LowNumCancerCell ): Label.append('SC') # significant control

    dic = {'sample':df['sample'],'replicate':df['replicate'],'imagens':listImagensNorm,'live_cancer':df['live_cancer'],'label':Label}
    new_df = pd.DataFrame(dic).reset_index(drop=True)
    return new_df

def Classifier_Patients(df_lastFrame):
    # Check how many samples are unique identifiable
    classeSample = ['NC','MC','SC','NC+MC','NC+SC','MC+SC','NC+MC+SC']
    sampleID = []
    labelSample = []
    meanImagens = []
    for id in df_lastFrame['sample'].unique():
        temp_df = df_lastFrame.loc[df_lastFrame['sample'] == id]
        sampleID.append(id)
        meanImagens.append(temp_df['imagens'].mean())
        Num_NC = len(temp_df.loc[temp_df['label'] == 'NC'])
        Num_MC = len(temp_df.loc[temp_df['label'] == 'MC'])
        Num_SC = len(temp_df.loc[temp_df['label'] == 'SC'])
        if ( len(temp_df['label'].unique()) == 1 ):
            if (Num_NC != 0 ): labelSample.append(classeSample[0])
            if (Num_MC != 0 ): labelSample.append(classeSample[1])
            if (Num_SC != 0 ): labelSample.append(classeSample[2])
        if ( len(temp_df['label'].unique()) == 2 ):
            if ( (Num_NC != 0) & (Num_MC != 0) ): labelSample.append(classeSample[3])
            if ( (Num_NC != 0) & (Num_SC != 0) ): labelSample.append(classeSample[4])
            if ( (Num_MC != 0) & (Num_SC != 0) ): labelSample.append(classeSample[5])
        if ( len(temp_df['label'].unique()) == 3 ):
            labelSample.append(classeSample[6])
    data_df_7classes = pd.DataFrame({'sample': sampleID,'imagens_mean':meanImagens,'label':labelSample})
    data_df_3classes = data_df_7classes.replace({'label': {'NC+MC':'Mixed', 'NC+SC':'Mixed', 'MC+SC':'Mixed', 'NC+MC+SC':'Mixed'}})
    return data_df_7classes, data_df_3classes

def Normalize_Parameters(df_input, Pat_df_3classes):
    df_input_min_max_scaled = pd.concat([df_input, Pat_df_3classes['label']], axis=1)
    # Normalize the parameters
    for column in NameParameters:
      MinValue = df_input_min_max_scaled[column].min()
      MaxValue = df_input_min_max_scaled[column].max()
      Delta_Value = MaxValue - MinValue
      df_input_min_max_scaled[column] = (df_input_min_max_scaled[column] - MinValue) / Delta_Value
    return df_input_min_max_scaled

def T_test(df_input_min_max_scaled):
    df_SC_NC_min_max_scaled = df_input_min_max_scaled.loc[(df_input_min_max_scaled["label"] == 'NC') | (df_input_min_max_scaled["label"] == 'SC')]
    Parameters_data = df_SC_NC_min_max_scaled
    T_test = []
    for parName in NameParameters:
        T_test.append( ttest_ind(Parameters_data.loc[Parameters_data['label'] == 'NC'][parName], Parameters_data.loc[Parameters_data['label'] == 'SC'][parName] ) )
    dict_tTest = {NameParameters[i]: T_test[i].statistic for i in range(len(NameParameters))}
    dict_pValue = {NameParameters[i]: T_test[i].pvalue for i in range(len(NameParameters))}
    return pd.DataFrame([dict_tTest, dict_pValue],index=['t-stat', 'p-value'])
    
def MannWhitneyU_test(df_input_min_max_scaled):
    df_SC_NC_min_max_scaled = df_input_min_max_scaled.loc[(df_input_min_max_scaled["label"] == 'NC') | (df_input_min_max_scaled["label"] == 'SC')]
    Parameters_data = df_SC_NC_min_max_scaled
    MWU_test = []
    for parName in NameParameters:
        MWU_test.append( mannwhitneyu(Parameters_data.loc[Parameters_data['label'] == 'NC'][parName], Parameters_data.loc[Parameters_data['label'] == 'SC'][parName], method="asymptotic" ) )
    dict_MWUTest = {NameParameters[i]: MWU_test[i].statistic for i in range(len(NameParameters))}
    dict_pValue = {NameParameters[i]: MWU_test[i].pvalue for i in range(len(NameParameters))}
    return pd.DataFrame([dict_MWUTest, dict_pValue],index=['MWU-stat', 'p-value'])

def Fitter_dist(df_input_min_max_scaled):
    dic_Fit = {}
    fileName = 'FitterPars.pickle'
    if not os.path.exists(fileName):
        df_SC_NC_min_max_scaled = df_input_min_max_scaled.loc[(df_input_min_max_scaled["label"] == 'NC') | (df_input_min_max_scaled["label"] == 'SC')]
        data_plot = pd.melt(df_SC_NC_min_max_scaled, id_vars=['sample','label'], value_vars=NameParameters)
        for parName in ['DC_max_recruitment_rate','macrophage_max_recruitment_rate']:
            f_SC = Fitter(data_plot.loc[(data_plot['label'] == 'SC') & (data_plot['variable'] == parName)]['value'].to_numpy()); f_SC.fit()
            f_NC = Fitter(data_plot.loc[(data_plot['label'] == 'NC') & (data_plot['variable'] == parName)]['value'].to_numpy()); f_NC.fit()
            dic_Fit[parName+'SC'] = f_SC
            dic_Fit[parName+'NC'] = f_NC
        with open(fileName, 'wb') as handle:
            pickle.dump(dic_Fit, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        dic_Fit = loading_pkl_file(fileName)
    return dic_Fit

def PatientsAnalysis():
    df = Loading_AditionalDataset()
    Mean_liveCancer_NC = []; Std_liveCancer_NC = []; ParMac_NC = []; ParDC_NC = []; Patient_NC = []; PatientLabel_NC = []
    Mean_liveCancer_MC = []; Std_liveCancer_MC = []; ParMac_MC = []; ParDC_MC = []; Patient_MC = []; PatientLabel_MC = []
    for sample in df['sample'].unique():
        if ( df.loc[(df['sample'] == sample) & (df['replicate'] == 0)]['sample_ref'].to_list() in  NC_Patients_samples): # (No control patients)
            Mean_liveCancer_NC.append(np.array(df.loc[df['sample'] == sample]['cancer_cell_day10']).mean())
            Std_liveCancer_NC.append(np.array(df.loc[df['sample'] == sample]['cancer_cell_day10']).std())
            ParMac_NC.append(df.loc[(df['sample'] == sample) & (df['replicate'] == 0)]['macrophage_max_recruitment_rate'])
            ParDC_NC.append(df.loc[(df['sample'] == sample) & (df['replicate'] == 0)]['DC_max_recruitment_rate'])
            Patient_NC.append(df.loc[(df['sample'] == sample) & (df['replicate'] == 0)]['sample_ref'])
            labels = df.loc[df['sample'] == sample]['label'].unique()
            if (len(labels) > 1) : PatientLabel_NC.append('Mixed')
            elif(labels[0] == 'NC'): PatientLabel_NC.append('NC')
            elif(labels[0] == 'SC'): PatientLabel_NC.append('SC')
            elif(labels[0] == 'MC'): PatientLabel_NC.append('MC')
        else: # (Marginal control patients)
            Mean_liveCancer_MC.append(np.array(df.loc[df['sample'] == sample]['cancer_cell_day10']).mean())
            Std_liveCancer_MC.append(np.array(df.loc[df['sample'] == sample]['cancer_cell_day10']).std())
            ParMac_MC.append(df.loc[(df['sample'] == sample) & (df['replicate'] == 0)]['macrophage_max_recruitment_rate'])
            ParDC_MC.append(df.loc[(df['sample'] == sample) & (df['replicate'] == 0)]['DC_max_recruitment_rate'])
            Patient_MC.append(df.loc[(df['sample'] == sample) & (df['replicate'] == 0)]['sample_ref'])
            labels = df.loc[df['sample'] == sample]['label'].unique()
            if (len(labels) > 1) : PatientLabel_MC.append('Mixed')
            elif(labels[0] == 'NC'): PatientLabel_MC.append('NC')
            elif(labels[0] == 'SC'): PatientLabel_MC.append('SC')
            elif(labels[0] == 'MC'): PatientLabel_MC.append('MC')

    labelInt = np.array(PatientLabel_NC+PatientLabel_MC) # concatenate the NC and MC patients
    labelInt[labelInt == 'NC'] = 0
    labelInt[labelInt == 'MC'] = 1; labelInt[labelInt == 'Mixed'] = 1 # combine Mixed and MC
    labelInt[labelInt == 'SC'] = 2
    return pd.DataFrame({'patientID': np.array(Patient_NC+Patient_MC).flatten(), 'Mac_recruit': np.array(ParMac_NC+ParMac_MC).flatten(), 'DC_recruit': np.array(ParDC_NC+ParDC_MC).flatten(), 'Mean_liveCancer': np.array(Mean_liveCancer_NC+Mean_liveCancer_MC), 'Std_liveCancer': np.array(Std_liveCancer_NC+Std_liveCancer_MC), 'label': labelInt})
