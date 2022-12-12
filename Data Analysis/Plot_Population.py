import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from DataSet import *

def plotSamplesLN(DataFrame, sampleID, ax):
  df = DataFrame.loc[DataFrame['sample'] == sampleID]
  df.loc[:]['simulation_time'] = df['simulation_time']/1440.0 # convert min to day
  dic1 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['DC'], 'Cells':len(df)*['DC']}
  dic2 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['TC'], 'Cells':len(df)*['TC']}
  dic3 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['TH1'], 'Cells':len(df)*['TH1']}
  dic4 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['TH2'], 'Cells':len(df)*['TH2']}
  dic5 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['TCt'], 'Cells':len(df)*['TCt']}
  dic6 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['Tht'], 'Cells':len(df)*['Tht']}
  new_df =pd.concat([pd.DataFrame(dic1), pd.DataFrame(dic2), pd.DataFrame(dic3), pd.DataFrame(dic4), pd.DataFrame(dic5), pd.DataFrame(dic6)],ignore_index=True)
  colorblindpallete = ['#016ca4','#cd5000','#ffbb7b','#595959','#ababab','#a1c8f0']
  lp = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', palette=colorblindpallete, data=new_df,ax=ax) # mean and 95% confidence interval
  lp.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
        ncol=3, fancybox=True, shadow=True)
  ax.set(yscale='symlog')

def plotSamplesMicroEnv(DataFrame, sampleID, ax1, ax2):
  df = DataFrame.loc[DataFrame['sample'] == sampleID]
  df.loc[:]['simulation_time'] = df['simulation_time']/1440.0 # convert min to day
  dic1 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['live_lung'], 'Cells':len(df)*['live_lung']}
  dic2 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['live_cancer'], 'Cells':len(df)*['live_cancer']}
  dic3 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['inac_DC'], 'Cells':len(df)*['inac_DC']}
  dic4 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['act_DC'], 'Cells':len(df)*['act_DC']}
  dic5 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['inac_Mac'], 'Cells':len(df)*['inac_Mac']}
  dic6 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['act_Mac'], 'Cells':len(df)*['act_Mac']}
  dic7 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['exas_Mac'], 'Cells':len(df)*['exas_Mac']}
  dic8 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['hyper_MAC'], 'Cells':len(df)*['hyper_MAC']}
  dic9 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['live_CD4'], 'Cells':len(df)*['live_CD4']}
  dic10 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['live_CD8'], 'Cells':len(df)*['live_CD8']}
  new_df1 =pd.concat([pd.DataFrame(dic1), pd.DataFrame(dic2)],ignore_index=True)
  CD8_mean = np.array(df.loc[df['simulation_time'] == 10.0]['live_CD8']).mean()
  ActMac_mean = np.array(df.loc[df['simulation_time'] == 10.0]['live_CD8']).mean()
  # Concatenate dictionaries
  new_df2 = pd.concat([pd.DataFrame(dic3), pd.DataFrame(dic4), pd.DataFrame(dic5), pd.DataFrame(dic6),pd.DataFrame(dic7), pd.DataFrame(dic8), pd.DataFrame(dic9), pd.DataFrame(dic10)],ignore_index=True)
  colours = {'live_lung': 'b', 'live_cancer': 'y', 'inac_DC': '#810F7C', 'act_DC': '#ff1493', 'inac_Mac': '#238B45', 'act_Mac': '#C0FF00', 'exas_Mac': '#A8DD76', 'hyper_MAC': '#A8DDB5', 'live_CD4': 'orange','live_CD8': 'red'}

  CD8_mean_threshold = 0
  ActMac_mean_threshold = 0
  if ( CD8_mean >= CD8_mean_threshold and ActMac_mean >= ActMac_mean_threshold):
    new_df2_1 = pd.concat([pd.DataFrame(dic3), pd.DataFrame(dic4), pd.DataFrame(dic5),pd.DataFrame(dic7), pd.DataFrame(dic8), pd.DataFrame(dic9)],ignore_index=True)
    new_df2_2 = pd.concat([pd.DataFrame(dic6), pd.DataFrame(dic10)],ignore_index=True)
    lp2_1 = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', data=new_df2_1,palette=[colours['inac_DC'],colours['act_DC'],colours['inac_Mac'],colours['exas_Mac'],colours['hyper_MAC'],colours['live_CD4']],ax=ax2) # mean and 95% confidence interval
    ax1_2 = ax2.twinx()
    ax1_2.set_ylabel('Number of Cells',color='blue')
    ax1_2.tick_params(axis='y',color='blue', labelcolor='blue')
    lp2_2 = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', data=new_df2_2,palette=[colours['act_Mac'],colours['live_CD8']],ax=ax1_2)
    handles, labels = [(a + b) for a, b in zip(ax2.get_legend_handles_labels(), ax1_2.get_legend_handles_labels())]
    leg = ax2.legend(handles, labels,loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=4, fancybox=True, shadow=True)
    ax1_2.get_legend().remove()
    for text in leg.get_texts():
        if (text._text == 'act_Mac' or text._text == 'live_CD8'): text.set_color("blue")
  elif ( CD8_mean > CD8_mean_threshold ):
    new_df2_1 = pd.concat([pd.DataFrame(dic3), pd.DataFrame(dic4), pd.DataFrame(dic5),pd.DataFrame(dic6), pd.DataFrame(dic7), pd.DataFrame(dic8), pd.DataFrame(dic9), ],ignore_index=True)
    new_df2_2 = pd.concat([pd.DataFrame(dic10)],ignore_index=True)
    lp2_1 = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', data=new_df2_1,palette=[colours['inac_DC'],colours['act_DC'],colours['inac_Mac'],colours['act_Mac'],colours['exas_Mac'],colours['hyper_MAC'],colours['live_CD4']],ax=ax2) # mean and 95% confidence interval
    ax1_2 = ax2.twinx()
    ax1_2.set_ylabel('Number of Cells',color='blue')
    ax1_2.tick_params(axis='y',color='blue', labelcolor='blue')
    lp2_2 = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', data=new_df2_2,palette=[colours['live_CD8']],ax=ax1_2)
    handles, labels = [(a + b) for a, b in zip(ax2.get_legend_handles_labels(), ax1_2.get_legend_handles_labels())]
    leg = ax2.legend(handles, labels,loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=4, fancybox=True, shadow=True)
    ax1_2.get_legend().remove()
    for text in leg.get_texts():
        if (text._text == 'live_CD8'): text.set_color("blue")
  elif ( ActMac_mean > ActMac_mean_threshold ):
    new_df2_1 = pd.concat([pd.DataFrame(dic3), pd.DataFrame(dic4), pd.DataFrame(dic5), pd.DataFrame(dic7), pd.DataFrame(dic8), pd.DataFrame(dic9), pd.DataFrame(dic10)],ignore_index=True)
    new_df2_2 = pd.concat([pd.DataFrame(dic6)],ignore_index=True)
    lp2_1 = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', data=new_df2_1,palette=[colours['inac_DC'],colours['act_DC'],colours['inac_Mac'],colours['exas_Mac'],colours['hyper_MAC'],colours['live_CD4'],colours['live_CD8']],ax=ax2) # mean and 95% confidence interval
    ax1_2 = ax2.twinx()
    ax1_2.set_ylabel('Number of Cells',color='blue')
    ax1_2.tick_params(axis='y',color='blue', labelcolor='blue')
    lp2_2 = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', data=new_df2_2,palette=[colours['act_Mac']],ax=ax1_2)
    handles, labels = [(a + b) for a, b in zip(ax2.get_legend_handles_labels(), ax1_2.get_legend_handles_labels())]
    leg = ax2.legend(handles, labels,loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=4, fancybox=True, shadow=True)
    ax1_2.get_legend().remove()
    for text in leg.get_texts():
        if (text._text == 'act_Mac'): text.set_color("blue")
  else:
    lp2 = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', data=new_df2,palette=[colours['inac_DC'],colours['act_DC'],colours['inac_Mac'],colours['act_Mac'],colours['exas_Mac'],colours['hyper_MAC'],colours['live_CD4'],colours['live_CD8']],ax=ax2) # mean and 95% confidence interval
    lp2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
          ncol=4, fancybox=True, shadow=True)
  lp1 = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', data=new_df1,palette=[colours['live_lung'],colours['live_cancer']],ax=ax1) # mean and 95% confidence interval
  lp1.axhline(100,ls='-.',c='gray',label='Threshold for NC')
  lp1.axhline(1500,ls='--',c='gray',label='Threshold for SC')
  lp1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
        ncol=2, fancybox=True, shadow=True)

def plot_PopCells(df_output, sampleID): # Patients examples: NC=2, MC = 93, SC=5, NC+MC=9, NC+SC=7, MC+SC=106, -- (CC and MC and SC)= 6
    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(16, 5))
    plotSamplesLN(df_output, sampleID,ax1)
    plotSamplesMicroEnv(df_output, sampleID,ax2,ax3)
    # set the spacing between subplots
    fig.tight_layout()
    plt.show()
