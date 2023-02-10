import matplotlib.pyplot as plt
import numpy as np
import glob

def annotate_axes(ax, text, fontsize=18):
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", fontsize=fontsize, color="darkgrey")

def PlotFrame(popcells, imagecell, frameID, FigName=None):
  # ax0: Tumor microenvironment
  # ax1: Legend
  # ax2: Melanoma and healthy tissue population
  # ax3: Lymph Node population
  # ax4: Immune cells population
  fig = plt.figure(layout=None, facecolor='0.9',figsize=(12,12))
  gs = fig.add_gridspec(nrows=3, ncols=3, top=0.94, bottom=0.1, left=0.05, right=0.95, hspace=0.3, wspace=0.3)
  ax0 = fig.add_subplot(gs[:-1, :-1])
  cells = plt.imread(imagecell)
  ax0.imshow(cells)
  ax0.axes.xaxis.set_visible(False)
  ax0.axes.yaxis.set_visible(False)
  ax1 = fig.add_subplot(gs[:-1, -1])
  legend = plt.imread('../beta/legend.png')
  ax1.imshow(legend)
  ax1.axes.xaxis.set_visible(False)
  ax1.axes.yaxis.set_visible(False)
  # Population plots
  colours = {'live_lung': 'b', 'live_cancer': 'y', 'inac_DC': '#810F7C', 'act_DC': '#ff1493', 'inac_Mac': '#238B45', 'act_Mac': '#C0FF00', 'exas_Mac': '#A8DD76', 'hyper_MAC': '#A8DDB5', 'live_CD4': 'orange','live_CD8': 'red', 'LN_DC':'#016ca4', 'LN_TC':'#cd5000', 'LN_TH1':'#ffbb7b', 'LN_TH2':'#595959', 'LN_TCt':'#ababab', 'LN_Tht':'#a1c8f0'}
  # Melanoma and Lung
  ax2 = fig.add_subplot(gs[-1, 0])
  ax2.plot(popcells[:,0]/1440.0, popcells[:,7], color=colours['live_lung'], label='live_lung')
  ax2.plot(popcells[:,0]/1440.0, popcells[:,14], color=colours['live_cancer'], label='live_cancer')
  ax2.scatter(popcells[frameID,0]/1440.0, popcells[frameID,7], color=colours['live_lung'])
  ax2.scatter(popcells[frameID,0]/1440.0, popcells[frameID,14], color=colours['live_cancer'])
  ax2.set(xlabel='Time (days)', ylabel='Number of cells')
  ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=2, fancybox=True, shadow=True, prop={'size': 8})
  ax2.ticklabel_format(axis='y', style='sci', scilimits=(-1,1), useMathText=True)
  # Lymph Node
  ax3 = fig.add_subplot(gs[-1, 1])
  ax3.set(yscale='symlog')
  ax3.plot(popcells[:,0]/1440.0, popcells[:,1], color=colours['LN_DC'], label='DC')
  ax3.scatter(popcells[frameID,0]/1440.0, popcells[frameID,1], color=colours['LN_DC'])
  ax3.plot(popcells[:,0]/1440.0, popcells[:,2], color=colours['LN_TC'], label='TC')
  ax3.scatter(popcells[frameID,0]/1440.0, popcells[frameID,2], color=colours['LN_TC'])
  ax3.plot(popcells[:,0]/1440.0, popcells[:,3], color=colours['LN_TH1'], label='TH1')
  ax3.scatter(popcells[frameID,0]/1440.0, popcells[frameID,3], color=colours['LN_TH1'])
  ax3.plot(popcells[:,0]/1440.0, popcells[:,4], color=colours['LN_TH2'], label='TH2')
  ax3.scatter(popcells[frameID,0]/1440.0, popcells[frameID,4], color=colours['LN_TH2'])
  ax3.plot(popcells[:,0]/1440.0, popcells[:,5], color=colours['LN_TCt'], label='TCt')
  ax3.scatter(popcells[frameID,0]/1440.0, popcells[frameID,5], color=colours['LN_TCt'])
  ax3.plot(popcells[:,0]/1440.0, popcells[:,6], color=colours['LN_Tht'], label='Tht')
  ax3.scatter(popcells[frameID,0]/1440.0, popcells[frameID,6], color=colours['LN_Tht'])
  ax3.set(xlabel='Time (days)', ylabel='Number of cells (log scale)')
  ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=3, fancybox=True, shadow=True, prop={'size': 8})
  # Immune cells
  ax4 = fig.add_subplot(gs[-1, 2])
  ax4_2 = ax4.twinx()
  ax4_2.set_ylabel('Number of Cells',color='blue')
  ax4_2.ticklabel_format(axis='y', style='sci', scilimits=(-1,1), useMathText=True)
  ax4_2.tick_params(axis='y',color='blue', labelcolor='blue')
  ax4.plot(popcells[:,0]/1440.0, popcells[:,15], color=colours['inac_DC'],label='inac_DC')
  ax4.scatter(popcells[frameID,0]/1440.0, popcells[frameID,15], color=colours['inac_DC'])
  ax4.plot(popcells[:,0]/1440.0, popcells[:,20], color=colours['act_DC'],label='act_DC')
  ax4.scatter(popcells[frameID,0]/1440.0, popcells[frameID,20], color=colours['act_DC'])
  ax4.plot(popcells[:,0]/1440.0, popcells[:,16], color=colours['inac_Mac'],label='inac_Mac')
  ax4.scatter(popcells[frameID,0]/1440.0, popcells[frameID,16], color=colours['inac_Mac'])
  ax4_2.plot(popcells[:,0]/1440.0, popcells[:,21], color=colours['act_Mac'],label='act_Mac')
  ax4_2.scatter(popcells[frameID,0]/1440.0, popcells[frameID,21], color=colours['act_Mac'])
  ax4.plot(popcells[:,0]/1440.0, popcells[:,17], color=colours['exas_Mac'],label='exas_Mac')
  ax4.scatter(popcells[frameID,0]/1440.0, popcells[frameID,17], color=colours['exas_Mac'])
  ax4.plot(popcells[:,0]/1440.0, popcells[:,18], color=colours['hyper_MAC'],label='hyper_MAC')
  ax4.scatter(popcells[frameID,0]/1440.0, popcells[frameID,18], color=colours['hyper_MAC'])
  ax4.plot(popcells[:,0]/1440.0, popcells[:,19], color=colours['live_CD4'],label='live_CD4')
  ax4.scatter(popcells[frameID,0]/1440.0, popcells[frameID,19], color=colours['live_CD4'])
  ax4_2.plot(popcells[:,0]/1440.0, popcells[:,22], color=colours['live_CD8'],label='live_CD8')
  ax4_2.scatter(popcells[frameID,0]/1440.0, popcells[frameID,22], color=colours['live_CD8'])
  ax4.set(xlabel='Time (days)', ylabel='Number of cells')
  ax4.ticklabel_format(axis='y', style='sci', scilimits=(-1,1), useMathText=True)
  # Concatenate legends
  handles, labels = [(a + b) for a, b in zip(ax4.get_legend_handles_labels(), ax4_2.get_legend_handles_labels())]
  leg = ax4.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=4, fancybox=True, shadow=True, prop={'size': 7})
  for text in leg.get_texts():
      if (text._text == 'act_Mac' or text._text == 'live_CD8'): text.set_color("blue")
  if (FigName):
      plt.savefig(FigName, bbox_inches = 'tight', pad_inches = 0)
      plt.close()
  else: plt.show()


if __name__ == '__main__':
    folderPatient = ['outputS000002','outputS000005','outputS000093']
    for folder in folderPatient:
        popcells = np.loadtxt('../Results/'+folder+'/NumCells.dat')
        dir_path = '../Results/'+folder+'/snapshot*.jpg'
        imagecell = glob.glob(dir_path)
        for frameID in range(len(imagecell)):
            PlotFrame(popcells,imagecell[frameID],frameID,f'../Results/%s/Video/frame_%.4d.jpg'%(folder,frameID))
