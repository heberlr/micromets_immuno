import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib_venn import venn3_unweighted
import seaborn as sns
import numpy as np
import pandas as pd
from DataSet import *

colours = {'NC':'#CA8151','MC':'#FFE1B3','SC':'#307CA1','NC+MC':'#999999','NC+SC':'#797979','MC+SC':'#898989','NC+MC+SC':'#696969','Mixed':'#595959', 'Mixed/MC': '#494949'}

def Plot_Patients(Sim_df_lastFrame,Pat_df_7classes,Pat_df_3classes, Title=False, FigName=None):
  fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(16, 6))
  Sim_sizes = Sim_df_lastFrame['label'].value_counts()
  Pat_7classes_sizes =  Pat_df_7classes['label'].value_counts()
  Pat_3classes_sizes =  Pat_df_3classes['label'].value_counts()
  ax1.pie(Sim_sizes,labels=Sim_sizes.index,colors=[colours[key] for key in Sim_sizes.index], autopct='%1.1f%%', shadow=True, startangle=90, textprops={'fontsize': 14,'color':'k'})
  if (Title): ax1.set_title("All 100k simulations")
  Venn = venn3_unweighted(subsets=(Pat_7classes_sizes['NC'],Pat_7classes_sizes['MC'],Pat_7classes_sizes['NC+MC'],Pat_7classes_sizes['SC'],Pat_7classes_sizes['NC+SC'],Pat_7classes_sizes['MC+SC'],Pat_7classes_sizes['NC+MC+SC']), set_labels=('NC','MC','SC','NC+MC','NC+SC','MC+SC','NC+MC+SC'),set_colors=(colours['NC'],colours['MC'],colours['SC']),subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/10000.0 ):1.1%}" + ")", alpha = 1,ax=ax2)
  Venn.get_patch_by_id('110').set_color(colours['NC+MC'])
  Venn.get_patch_by_id('101').set_color(colours['NC+SC'])
  Venn.get_patch_by_id('011').set_color(colours['MC+SC'])
  Venn.get_patch_by_id('111').set_color(colours['NC+MC+SC'])
  for text in Venn.set_labels:
     text.set_fontsize(14)
  for text in Venn.subset_labels:
     text.set_fontsize(12)
  if (Title): ax2.set_title("10k patients - 7 sets")
  ax3.pie(Pat_3classes_sizes,labels=Pat_3classes_sizes.index,colors=[colours[key] for key in Pat_3classes_sizes.index], autopct='%1.1f%%', shadow=True, startangle=90, textprops={'fontsize': 14,'color':'k'})
  if (Title): ax3.set_title("10k patients - 3 sets")
  if (FigName): plt.savefig(FigName)
  else: plt.show()

def Plot_Parameters(df_input_min_max_scaled, FigName=None):
  df_SC_NC_min_max_scaled = df_input_min_max_scaled.loc[(df_input_min_max_scaled["label"] == 'NC') | (df_input_min_max_scaled["label"] == 'SC')]
  data_plot = pd.melt(df_SC_NC_min_max_scaled, id_vars=['sample','label'], value_vars=NameParameters)
  plt.figure(figsize=(16, 6))
  ax = sns.boxplot(x="variable", y="value", hue="label", data=data_plot, palette=colours, saturation=1)
  ax.set_xticklabels([r"$r_{recruit}[MP]$",r"$\rho_{min}[MP]$",r"$\rho_{sat}[MP]$",r"$r_{recruit}[DC]$",r"$\rho_{min}[DC]$",r"$\rho_{sat}[DC]$",r"$r_{leave}$",r"$\delta_C$",r"$\kappa_{T}$",r"$\delta_{DM}$"],fontsize=14)
  ax.set(xlabel=None)
  ax.set_ylabel('Normalized parameters',fontsize=14)
  ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1),fontsize=14)
  if (FigName): plt.savefig(FigName)
  else: plt.show()

def Plot_PatientAnalysis(PlotPatientPoint = False, FigName=None):
    df_PatientVariation = PatientsAnalysis()
    fig = plt.figure(figsize=(12, 10))
    G = gridspec.GridSpec(4, 2, width_ratios=[1,0.05],height_ratios=[0.01,1,0.01,1])
    # Patient classified as marginal control
    ax_title_1 = plt.subplot(G[0, 0])
    ax_title_1.axis('off')
    ax_title_1.set_title('Ten virtual patients initially classified as marginal control (MC)')
    G0 = gridspec.GridSpecFromSubplotSpec(2, 5, subplot_spec=G[1,0])
    # Patient classified as no control
    ax_title_2 = plt.subplot(G[2, 0])
    ax_title_2.axis('off')
    ax_title_2.set_title('Ten virtual patients initially classified as no control (NC)')
    G1 = gridspec.GridSpecFromSubplotSpec(2, 5, subplot_spec=G[3,0])
    # Axes color bar
    axcb = plt.subplot(G[:, 1])
    NewColours = {'NC': colours['NC'], 'Mixed/MC': colours['Mixed/MC'], 'SC': colours['SC']}
    customPalette = sns.color_palette(NewColours.values(),len(NewColours))

    for ind, patientID in enumerate(MC_Patients_samples):
        ax = plt.subplot(G0[int(ind/5), ind%5])
        df = df_PatientVariation.loc[df_PatientVariation['patientID'] == patientID]
        df_wide = df.pivot_table( index='DC_recruit', columns='Mac_recruit', values='label')
        g1 = sns.heatmap(df_wide,cmap=customPalette,vmin=0,vmax=2,cbar=False,ax=ax,linewidths=.5)#,annot=True)
        g1.invert_yaxis()
        g1.set(xlabel = (''), ylabel = (''), xticks = ([]), yticks = ([]))
        g1.set_title('Patient '+str(patientID))#,color="#4783af")
        if ( ind%5 == 0): g1.set_yticks(ticks=[0.01,1.5,2.99],labels=['0','0.5','1'])
        if (PlotPatientPoint):
            MAC_rec_value = df_input.loc[df_input['sample'] == patientID]['macrophage_max_recruitment_rate'].values[0]
            DC_rec_value = df_input.loc[df_input['sample'] == patientID]['DC_max_recruitment_rate'].values[0]
            ax.scatter(3*MAC_rec_value/8e-9,3*DC_rec_value/4e-9, color=colours['MC'], edgecolors='k')

    for ind, patientID  in enumerate(NC_Patients_samples):
        ax = plt.subplot(G1[int((ind)/5), ind%5])
        df = df_PatientVariation.loc[df_PatientVariation['patientID'] == patientID]
        df_wide = df.pivot_table( index='DC_recruit', columns='Mac_recruit', values='label')
        if (ind == 9): g1 = sns.heatmap(df_wide,cmap=customPalette,vmin=0,vmax=2,ax=ax,linewidths=.5,cbar_ax=axcb)#,annot=True)
        else: g1 = sns.heatmap(df_wide,cmap=customPalette,vmin=0,vmax=2,cbar=False,ax=ax,linewidths=.5)#,annot=True)
        g1.invert_yaxis()
        g1.set(xlabel = (''), ylabel = (''), xticks = ([]), yticks = ([]))
        g1.set_title('Patient '+str(patientID))#,color="#ff7f7f")
        if ( int(ind/5) == 1): g1.set_xticks(ticks=[0.02,1.5,2.98],labels=['0','0.5','1'])
        if ( ind%5 == 0): g1.set_yticks(ticks=[0.01,1.5,2.99],labels=['0','0.5','1'])
        if (PlotPatientPoint):
            MAC_rec_value = df_input.loc[df_input['sample'] == patientID]['macrophage_max_recruitment_rate'].values[0]
            DC_rec_value = df_input.loc[df_input['sample'] == patientID]['DC_max_recruitment_rate'].values[0]
            ax.scatter(3*MAC_rec_value/8e-9,3*DC_rec_value/4e-9, color=colours['NC'], edgecolors='k')

    fig.supxlabel(r'$r_{recruit}[MP]$')
    fig.supylabel(r'$r_{recruit}[DC]$')
    # modify colorbar:
    colorbar = ax.collections[0].colorbar
    r = colorbar.vmax - colorbar.vmin
    NumPoints = len(NewColours.keys())
    colorbar.set_ticks([colorbar.vmin + r / NumPoints * (0.5 + i) for i in range(NumPoints)])
    colorbar.set_ticklabels(NewColours.keys())
    plt.tight_layout()
    if (FigName): plt.savefig(FigName)
    else: plt.show()

def Plot_UMAP(Sim_df_lastFrame, FigName=None):
    fig = plt.figure(figsize=(8, 7))
    df_UMAP, df = Loading_UMAP()
    df_UMAP['label'] = Sim_df_lastFrame['label']
    ax = sns.scatterplot(data=df_UMAP,x='UMAP1',y='UMAP2',hue='label',palette=[colours['NC'],colours['MC'],colours['SC']], hue_order=['NC','MC','SC'])
    ax.legend(markerscale=2.0, fontsize=14)
    ax.set_xlabel('UMAP 1', fontsize=14)
    ax.set_ylabel('UMAP 2', fontsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    if (FigName): plt.savefig(FigName)
    else: plt.show()

def Plot_UMAP_MeanPatients(Pat_df_7classes, FigName=None):
    fig = plt.figure(figsize=(8, 7))
    df,df_UMAP = Loading_UMAP()
    df_UMAP['label'] = Pat_df_7classes['label']
    ax = sns.scatterplot(data=df_UMAP,x='UMAP1',y='UMAP2',hue='label',palette=[colours['NC'],colours['MC'],colours['SC'],colours['NC+MC'],colours['NC+SC'],colours['MC+SC'],colours['NC+MC+SC']], hue_order=['NC','MC','SC','NC+MC','NC+SC','MC+SC','NC+MC+SC'])
    ax.legend(markerscale=2.0, fontsize=14)
    ax.set_xlabel('UMAP 1', fontsize=14)
    ax.set_ylabel('UMAP 2', fontsize=14)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    if (FigName): plt.savefig(FigName)
    else: plt.show()

if __name__ == '__main__':
    df_input, df_output = Loading_dataset()
    Sim_df_lastFrame = Classifier_Simulations(df_output)
    Pat_df_7classes, Pat_df_3classes = Classifier_Patients(Sim_df_lastFrame)
    Plot_UMAP(Sim_df_lastFrame,FigName='Figure5_D.svg') # Plot UMAP of the trajectories
    Plot_UMAP_MeanPatients(Pat_df_7classes,FigName='Figure5_B.svg') # Plot UMAP of the trajectory averages
    Plot_Patients(Sim_df_lastFrame,Pat_df_7classes,Pat_df_3classes, FigName='Figure5_AC.svg') # Patient statistics
    # Plot the quartiles of patient features
    df_input_min_max_scaled = Normalize_Parameters(df_input, Pat_df_3classes)
    Plot_Parameters(df_input_min_max_scaled, FigName='Figure4_A.svg')
    # Plot the analysis of 20 patients
    Plot_PatientAnalysis(FigName='Figure4_B.svg')
