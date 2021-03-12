from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import sys

class FormatScalarFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, fformat="%1.1f", offset=True, mathText=True):
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,
                                                        useMathText=mathText)
    def _set_format(self):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % self.format


if (len(sys.argv) != 5):
  print("Please provide 4 args: InitialTime LastTime, SavePNG (bool), and folder")
  sys.exit(1)
initial_index = int(sys.argv[1]);
last_index = int(sys.argv[2]);
SavePNG = int(sys.argv[3])
folder = sys.argv[4]
Lcell_size = 10;
Dcell_size = 5;
Ecell_size = 12;
skin_live_count = np.zeros( last_index+1 );
skin_live_mutated_count = np.zeros( last_index+1 );
skin_dead_count = np.zeros( last_index+1 );
macrophage_count = np.zeros( last_index+1 );
neuthophil_count = np.zeros( last_index+1 );
dendritic_count = np.zeros( last_index+1 );
CD8_count = np.zeros( last_index+1 );
CD4_count = np.zeros( last_index+1 );
times = np.zeros( last_index+1 );

figure, axes = plt.subplots(nrows=3, ncols=3,figsize=(10,8))

for n in range( initial_index,last_index+1 ):
  filename= "output%08i"%n+'.xml'
  filenameOut=folder+'/output'+"%08i"%n+'.png'
  mcds=pyMCDS(filename,folder)
  times[n]= mcds.get_time()

  cx = mcds.data['discrete_cells']['position_x'];
  cy = mcds.data['discrete_cells']['position_y'];
  cycle = mcds.data['discrete_cells']['cycle_model']
  current_phase = mcds.data['discrete_cells']['current_phase']
  current_phase = current_phase.astype(int)
  elapsed_time_in_phase = mcds.data['discrete_cells']['elapsed_time_in_phase']
  cell_type = mcds.data['discrete_cells']['cell_type']
  cell_type = cell_type.astype(int)

  DG = mcds.data['discrete_cells']['neoantigens_intracellular']

  skin_live_mutated = np.argwhere( (DG > 0.0) & (cycle < 100) & (cell_type==1) ).flatten()
  skin_live = np.argwhere( (DG == 0.0) & (cycle < 100) & (cell_type==1) ).flatten()
  skin_dead = np.argwhere( (cycle >= 100) & (cell_type==1) ).flatten()
  macrophage = np.argwhere( cell_type==4 ).flatten()
  neuthophil = np.argwhere( cell_type==5 ).flatten()
  dendritic = np.argwhere( cell_type==6 ).flatten()
  CD8 = np.argwhere( cell_type==3 ).flatten()
  CD4 = np.argwhere( cell_type==7 ).flatten()

  skin_live_count[n] = len(skin_live)
  skin_live_mutated_count[n] = len(skin_live_mutated)
  skin_dead_count[n] = len(skin_dead)
  macrophage_count[n] = len(macrophage)
  neuthophil_count[n] = len(neuthophil)
  dendritic_count[n] = len(dendritic)
  CD8_count[n] = len(CD8)
  CD4_count[n] = len(CD4)

  RadiusSize = 400.0

  # figure.suptitle( '#NC:'+str("%04i"%(skin_live_count[n]))+'  #MC:'+str("%04i"%(skin_live_mutated_count[n]))+ '--  #M:'+str("%04i"%(macrophage_count[n]))+ '  #N:'+str("%04i"%(neuthophil_count[n]))+ '  #DC:'+str("%04i"%(dendritic_count[n]))+'\n'+'  #CD8:'+str("%04i"%(CD8_count[n]))+ '  #CD4:'+str("%04i"%(CD4_count[n]))+ '  #D:'+str("%04i"%(skin_dead_count[n]))+'  Time:' +str("%8.2f"%(n)) + ' hours', size=14)
  #
  # plt.subplot(221)
  # plt.scatter( cx[skin_live],cy[skin_live],c='blue',s=Lcell_size,label='NC');
  # plt.scatter( cx[skin_live_mutated],cy[skin_live_mutated],c='yellow',s=Lcell_size,label='MC');
  # plt.scatter( cx[skin_dead],cy[skin_dead],c='black',s=Dcell_size,label='D', alpha=0.25);
  # plt.scatter( cx[macrophage],cy[macrophage],c='green',s=Ecell_size,label='M' );
  # plt.scatter( cx[neuthophil],cy[neuthophil],c='cyan',s=Ecell_size,label='N' );
  # plt.scatter( cx[dendritic],cy[dendritic],c='brown',s=Ecell_size,label='DC' );
  # plt.scatter( cx[CD8],cy[CD8],c='red',s=Ecell_size,label='CD8' );
  # plt.scatter( cx[CD4],cy[CD4],c='orange',s=Ecell_size,label='CD4' );
  # plt.xlim(-RadiusSize, RadiusSize)
  # plt.ylim(-RadiusSize, RadiusSize)
  # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
  # plt.subplots_adjust(left=0.08,right=0.93,bottom=0.06,top=0.89,wspace=0.26,hspace=0.26)
  #
  # plt.subplot(333)
  # DG = mcds.get_concentrations( 'neoantigens' );
  # X1,Y1 = mcds.get_2D_mesh();
  # if (DG.max() > 0):
  #   # v1 = np.linspace(0, DG.max(), 10, endpoint=True)
  #   v1 = np.linspace(0, 1.0, 10, endpoint=True)
  #   plt.contourf(X1,Y1,DG[:,:,0],levels=v1,cmap='binary');
  #   cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.1f"))
  # else:
  #   plt.contourf(X1,Y1,DG[:,:,0],cmap='binary');
  #   cbar = plt.colorbar(format=FormatScalarFormatter("%.1f"))
  # plt.xlim(-RadiusSize, RadiusSize)
  # plt.ylim(-RadiusSize, RadiusSize)
  # plt.title("Neoantigens")
  #
  # plt.subplot(336)
  # chemokine = mcds.get_concentrations( 'chemokine' );
  # X1,Y1 = mcds.get_2D_mesh();
  # if (chemokine.max() > 0):
  #   # v1 = np.linspace(0, chemokine.max(), 10, endpoint=True)
  #   v1 = np.linspace(0, 1.0, 10, endpoint=True)
  #   plt.contourf(X1,Y1,chemokine[:,:,0],levels=v1,cmap='binary');
  #   cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.1f"))
  # else:
  #   plt.contourf(X1,Y1,chemokine[:,:,0],cmap='binary');
  #   x = plt.colorbar(format=FormatScalarFormatter("%.1f"))
  # plt.xlim(-RadiusSize, RadiusSize)
  # plt.ylim(-RadiusSize, RadiusSize)
  # plt.title("Chemokine")
  #
  #
  # plt.subplot(337)
  #
  # PIC = mcds.get_concentrations( 'pro-inflammatory cytokine' );
  # if (PIC.max() > 0):
  #   # v1 = np.linspace(0, PIC.max(), 10, endpoint=True)
  #   v1 = np.linspace(0, 1.0, 10, endpoint=True)
  #   plt.contourf(X1,Y1,PIC[:,:,0],v1,cmap='binary');
  #   cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.1f"))
  # else:
  #   plt.contourf(X1,Y1,PIC[:,:,0],cmap='binary');
  #   cbar = plt.colorbar(format=FormatScalarFormatter("%.1f"))
  # plt.xlim(-RadiusSize, RadiusSize)
  # plt.ylim(-RadiusSize, RadiusSize)
  # plt.title("Pro-inf cytokine")
  #
  # plt.subplot(338)
  #
  # debris = mcds.get_concentrations( 'debris' );
  # if (debris.max() > 0):
  #   v1 = np.linspace(0, 1.0, 10, endpoint=True)
  #   # v1 = np.linspace(0, debris.max(), 10, endpoint=True)
  #   plt.contourf(X1,Y1,debris[:,:,0],v1,cmap='binary');
  #   cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.1f"))
  # else:
  #   plt.contourf(X1,Y1,debris[:,:,0],cmap='binary');
  #   cbar = plt.colorbar(format=FormatScalarFormatter("%.1f"))
  # plt.xlim(-RadiusSize, RadiusSize)
  # plt.ylim(-RadiusSize, RadiusSize)
  # plt.title("Debris")
  #
  #
  # plt.subplot(339)
  #
  # AIC = mcds.get_concentrations( 'anti-inflammatory cytokine' );
  # if (AIC.max() > 0):
  #   # v1 = np.linspace(0, AIC.max(), 10, endpoint=True)
  #   v1 = np.linspace(0, 1.0, 10, endpoint=True)
  #   plt.contourf(X1,Y1,AIC[:,:,0],v1,cmap='binary');
  #   cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.1f"))
  # else:
  #   plt.contourf(X1,Y1,AIC[:,:,0],cmap='binary');
  #   cbar = plt.colorbar(format=FormatScalarFormatter("%.1f"))
  # plt.xlim(-RadiusSize, RadiusSize)
  # plt.ylim(-RadiusSize, RadiusSize)
  # plt.title("Anti-inf cytokine")
  #
  # if (SavePNG):
  #   figure.savefig(filenameOut)
  #   #plt.savefig(filenameOut)
  # else:
  #   plt.draw()
  #   plt.waitforbuttonpress(0) # this will wait for indefinite time
  #   #plt.pause(0.2)
  # plt.clf()

plt.figure(2)
# plt.plot(times,skin_live_count)
plt.plot(times/1440.0,skin_live_mutated_count, color='yellow',label='cancer')
plt.plot(times/1440.0,skin_dead_count, color='black',label='dead')
plt.plot(times/1440.0,macrophage_count, color='green',label='macrophage')
plt.plot(times/1440.0,neuthophil_count, color='cyan',label='neuthophil')
plt.plot(times/1440.0,dendritic_count, color='brown',label='dendritic')
plt.plot(times/1440.0,CD8_count, color='red',label='CD8')
#plt.plot(times,CD4_count, color='orange',label='CD4')
plt.xlabel("Time (days)")
plt.ylabel("Number of cell")
plt.legend(loc='upper left')
plt.show()
