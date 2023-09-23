from DataSet import *

class ScalarFormatterClass(ScalarFormatter):
   def _set_format(self):
      self.format = "%1.1f"

def plotSamplesLN(DataFrame, sampleID, ax, FontSize=None, replicateID = None, ymax=None):
  if (replicateID or replicateID == 0): df = DataFrame.loc[ (DataFrame['sample'] == sampleID) & (DataFrame['replicate'] == replicateID) ]
  else: df = DataFrame.loc[DataFrame['sample'] == sampleID]
  df.loc[:]['simulation_time'] = df['simulation_time']/1440.0 # convert min to day
  dic1 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['DC'], 'Cells':len(df)*['DC']}
  dic2 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['TC'], 'Cells':len(df)*['TC']}
  dic3 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['TH1'], 'Cells':len(df)*['TH1']}
  dic4 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['TH2'], 'Cells':len(df)*['TH2']}
  dic5 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['TCt'], 'Cells':len(df)*['TCt']}
  dic6 = {'Time (days)': df['simulation_time'], 'Number of Cells':df['Tht'], 'Cells':len(df)*['Tht']}
  Newlabels = {'DC': r'$D_M$', 'TC':r'$T_C$','TH1': r'$T_{H1}$', 'TH2':r'$T_{H2}$','TCt': r'$T_{Ct}$', 'Tht':r'$T_{Ht}$'}
  new_df = pd.concat([pd.DataFrame(dic1), pd.DataFrame(dic2), pd.DataFrame(dic3), pd.DataFrame(dic4), pd.DataFrame(dic5), pd.DataFrame(dic6)],ignore_index=True)
  colorblindpallete = ['#016ca4','#cd5000','#ffbb7b','#595959','#ababab','#a1c8f0']
  lp = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', palette=colorblindpallete, data=new_df,ax=ax) # mean and 95% confidence interval
  leg = lp.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
        ncol=3, fancybox=True, shadow=True)
  if (ymax): ax.set(yscale='symlog',ylim=(0,ymax))
  else: ax.set(yscale='symlog')
  if FontSize:
      for item in ([ax.title, ax.xaxis.label, ax.yaxis.label, ax.yaxis.offsetText] + ax.get_xticklabels() + ax.get_yticklabels()):
          item.set_fontsize(FontSize)
  for text in leg.get_texts():
      if (text._text in Newlabels.keys()): text._text = Newlabels[text._text] # change names of legend
      if FontSize: text.set_size(FontSize) # change the label size

def plotSamplesMicroEnv(DataFrame, sampleID, ax1, ax2, FontSize=None, replicateID = None):
  if (replicateID or replicateID == 0): df = DataFrame.loc[ (DataFrame['sample'] == sampleID) & (DataFrame['replicate'] == replicateID) ]
  else: df = DataFrame.loc[DataFrame['sample'] == sampleID]
  # print(df.info())
  df.loc[:,'simulation_time'] = df['simulation_time'].div(1440.0) # convert min to day
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
  Newlabels = {'live_lung': 'lung cells', 'live_cancer': 'cancer cells', 'inac_DC': 'iDC', 'act_DC': 'aDC', 'inac_Mac': 'M0', 'act_Mac': 'M1/M2', 'exas_Mac': 'eMP', 'hyper_MAC': 'hMP', 'live_CD4': 'CD4+','live_CD8': 'CD8+'}

  CD8_mean_threshold = 50
  ActMac_mean_threshold = 50
  ax1_2 = None
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
    leg = lp2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
          ncol=4, fancybox=True, shadow=True)
  
  yScalarFormatter = ScalarFormatterClass(useMathText=True)
  yScalarFormatter.set_powerlimits((-1,1))
  items = ([ax2.title, ax2.xaxis.label, ax2.yaxis.label, ax2.yaxis.offsetText]+ax2.get_xticklabels()+ ax2.get_yticklabels())
  if (ax1_2): 
    ax1_2.yaxis.set_major_formatter(yScalarFormatter)
    items += ([ax1_2.yaxis.label, ax1_2.yaxis.offsetText]+ ax1_2.get_yticklabels())
  if FontSize:
    for item in items:
        item.set_fontsize(FontSize)
    for text in leg.get_texts():
        text._text = Newlabels[text._text] # change names of legend
        if FontSize: text.set_size(FontSize) # change the label size


  lp1 = sns.lineplot(x='Time (days)', y='Number of Cells',hue='Cells',marker='o', data=new_df1,palette=[colours['live_lung'],colours['live_cancer']],ax=ax1) # mean and 95% confidence interval
  lp1.axhline(100,ls='-.',c='gray',label='Threshold for NC')
  lp1.axhline(1500,ls='--',c='gray',label='Threshold for SC')
  leg1 = lp1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
        ncol=2, fancybox=True, shadow=True)
  if FontSize:
      for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label, ax1.yaxis.offsetText] + ax1.get_xticklabels() + ax1.get_yticklabels()):
          item.set_fontsize(FontSize)
  for text in leg1.get_texts():
      if (text._text in Newlabels.keys()): text._text = Newlabels[text._text] # change names of legend
      if FontSize: text.set_size(FontSize) # change the label size
  ax1.yaxis.set_major_formatter(yScalarFormatter)

def plot_PopCells(df_output, sampleID,FontSize=None,FigName=None, Vertical=None, replicateID=None): # Patients examples: NC=2, MC = 93, SC=5, NC+MC=9, NC+SC=7, MC+SC=106, -- (CC and MC and SC)= 6
    if (Vertical): fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(4, 8))
    else: fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(16, 4))
    plotSamplesLN(df_output, sampleID,ax2,FontSize=FontSize,replicateID=replicateID)
    plotSamplesMicroEnv(df_output, sampleID,ax1,ax3,FontSize=FontSize,replicateID=replicateID)
    # set the spacing between subplots
    fig.tight_layout()
    if (FigName): plt.savefig(FigName)
    else: plt.show()

if __name__ == '__main__':
    df_input, df_output = Loading_dataset()
    # Plot of populations of patient 6
    # plot_PopCells(df_output, 6,FontSize=18, FigName='Figure6.svg')
    
    # Plot of analysis of trajectories of NC in UMAP
    # plot_PopCells(df_output, 2546,FontSize=12, replicateID=5, Vertical=True, FigName='Sample_2546_R05.svg')
    # plot_PopCells(df_output, 4225,FontSize=12, replicateID=8, Vertical=True,FigName='Sample_4225_R08.svg')
    # exit()
    # Plot multiples LymphNode dynamics of NC samples in UMAP
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(4, 9))
    plotSamplesLN(df_output, 2489,ax1,FontSize=12, replicateID=8, ymax=1e6)
    plotSamplesLN(df_output, 6555,ax2,FontSize=12, replicateID=3, ymax=1e6)
    plotSamplesLN(df_output, 7010,ax3,FontSize=12, replicateID=7, ymax=1e6)
    fig.tight_layout()
    plt.savefig('LNs_2489_6555_7010.svg')
    fig2, (ax2_1,ax2_2,ax2_3) = plt.subplots(3,1,figsize=(4, 9))
    plotSamplesLN(df_output, 5666,ax2_1,FontSize=12, replicateID=7, ymax=1e6)
    plotSamplesLN(df_output, 4967,ax2_2,FontSize=12, replicateID=0, ymax=1e6)
    plotSamplesLN(df_output, 8307,ax2_3,FontSize=12, replicateID=6, ymax=1e6)
    fig2.tight_layout()
    plt.savefig('LNs_5666_4967_8307.svg')
    # plt.show()
