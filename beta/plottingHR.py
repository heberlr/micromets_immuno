from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import seaborn as sns
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

class NumberOfCells:
    def __init__(self, sizeV):
        self.sizeV = sizeV
        self.times = np.zeros( last_index+1 );
        self.lung = np.zeros( last_index+1 );
        self.melanoma = np.zeros( last_index+1 );
        self.dead = np.zeros( last_index+1 );
        self.macrophage = np.zeros( last_index+1 )
        self.dendritic = np.zeros( last_index+1 );
        self.CD8 = np.zeros( last_index+1 );
        self.CD4 = np.zeros( last_index+1 );
    def __add__(self, other):
        Sum = NumberOfCells(self.sizeV)
        Sum.times = self.times
        Sum.lung = self.lung + other.lung
        Sum.melanoma = self.melanoma + other.melanoma
        Sum.dead = self.dead + other.dead
        Sum.macrophage = self.macrophage + other.macrophage
        Sum.dendritic = self.dendritic + other.dendritic
        Sum.CD8 = self.CD8 + other.CD8
        Sum.CD4 = self.CD4 + other.CD4
        return Sum
    def __sub__(self, other):
        Sum = NumberOfCells(self.sizeV)
        Sum.times = self.times
        Sum.lung = self.lung - other.lung
        Sum.melanoma = self.melanoma - other.melanoma
        Sum.dead = self.dead - other.dead
        Sum.macrophage = self.macrophage - other.macrophage
        Sum.dendritic = self.dendritic - other.dendritic
        Sum.CD8 = self.CD8 - other.CD8
        Sum.CD4 = self.CD4 - other.CD4
        return Sum
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)
    def __truediv__(self, other):
        if not isinstance(other, (int, float)):
            return NotImplemented
        self.lung = self.lung / other
        self.melanoma = self.melanoma / other
        self.dead = self.dead / other
        self.macrophage = self.macrophage / other
        self.dendritic = self.dendritic / other
        self.CD8 = self.CD8 / other
        self.CD4 = self.CD4 / other
        return self
    def __pow__(self, other):
        if not isinstance(other, (int, float)):
            return NotImplemented
        self.lung = self.lung ** other
        self.melanoma = self.melanoma ** other
        self.dead = self.dead ** other
        self.macrophage = self.macrophage ** other
        self.dendritic = self.dendritic ** other
        self.CD8 = self.CD8 ** other
        self.CD4 = self.CD4 ** other
        return self
    def __abs__(self):
        self.lung = abs(self.lung)
        self.melanoma = abs(self.melanoma)
        self.dead = abs(self.dead)
        self.macrophage = abs(self.macrophage)
        self.dendritic = abs(self.dendritic)
        self.CD8 = abs(self.CD8)
        self.CD4 = abs(self.CD4)
        return self

def Std_class(mean_QOI, list_QOI):
    std_QOI = NumberOfCells(mean_QOI.sizeV)
    for i in range(len(list_QOI)):
        std_QOI = std_QOI + (list_QOI[i] - mean_QOI)**2
    std_QOI = ( std_QOI / len(list_QOI) )**0.5
    return std_QOI

def Replicates_plot(initial_index,last_index,folder,Nsamples):
    list_QOI = [NumberOfCells(last_index-initial_index+1) for j in range(Nsamples)]
    for i in range(Nsamples):
        for n in range( initial_index,last_index+1 ):
            folderTemp = folder + '_sample'+"%i"%i
            filenameOut=folderTemp+'/output'+"%08i"%n+'.png'
            filename= "output%08i"%n+'.xml'
            mcds=pyMCDS(filename,folderTemp)
            list_QOI[i].times[n]= mcds.get_time()
            Read_QOIs(mcds,n,list_QOI[i])

    mean_QOI = sum(list_QOI) / len(list_QOI)
    std_QOI = Std_class(mean_QOI,list_QOI)

    plotCurves(mean_QOI,std_QOI)

def Neoantigens_plot(initial_index,last_index,folder):
    times = np.zeros( last_index-initial_index+1 );
    clonal_count = np.zeros( last_index-initial_index+1 );
    subclonal_count = np.zeros( last_index-initial_index+1 );
    shared_count = np.zeros( last_index-initial_index+1 );
    for n in range( initial_index,last_index+1 ):
        filenameOut=folder+'/output'+"%08i"%n+'.png'
        filename= "output%08i"%n+'.xml'
        mcds=pyMCDS(filename,folder)
        times[n]= mcds.get_time()
        TypeNeoantigen = mcds.data['discrete_cells']['neoantigen_type']
        cell_type = mcds.data['discrete_cells']['cell_type']
        cycle = mcds.data['discrete_cells']['cycle_model']
        Melanoma = np.argwhere( (cell_type==2)  & (cycle < 100) ).flatten()
        ClonalNeoantigen = np.argwhere( (cell_type==2)  & (cycle < 100) & (TypeNeoantigen == 0) ).flatten()
        SubClonalNeoantigen = np.argwhere( (cell_type==2) & (cycle < 100) & (TypeNeoantigen == 2) ).flatten()
        SharedNeoantigen = np.argwhere( (cell_type==2) & (cycle < 100) & (TypeNeoantigen == 1) ).flatten()
        clonal_count[n] = len(ClonalNeoantigen)/len(Melanoma)
        subclonal_count[n] = len(SubClonalNeoantigen)/len(Melanoma)
        shared_count[n] = len(SharedNeoantigen)/len(Melanoma)
    plt.plot(times/1440.0,clonal_count, 'o', color='blue', label='Clonal')
    plt.plot(times/1440.0,subclonal_count, 'o', color='red', label='Subclonal')
    plt.plot(times/1440.0,shared_count, 'o', color='yellow', label='Shared')
    plt.xlabel("Time (days)")
    plt.ylabel("Neoantigens fraction")
    plt.legend()
    plt.show()

def Read_QOIs(mcds,n,QOIs):
    cycle = mcds.data['discrete_cells']['cycle_model']
    cell_type = mcds.data['discrete_cells']['cell_type']
    cell_type = cell_type.astype(int)

    melanoma = np.argwhere( (cell_type==2) & (cycle < 100) ).flatten()
    lung = np.argwhere( (cell_type==1) & (cycle < 100) ).flatten()
    dead = np.argwhere( (cycle >= 100) & ((cell_type==1) | (cell_type==2)) ).flatten()
    macrophage = np.argwhere( cell_type==4 ).flatten()
    dendritic = np.argwhere( cell_type==6 ).flatten()
    CD8 = np.argwhere( cell_type==3 ).flatten()
    CD4 = np.argwhere( cell_type==7 ).flatten()

    QOIs.lung[n] = len(lung)
    QOIs.melanoma[n] = len(melanoma)
    QOIs.dead[n] = len(dead)
    QOIs.macrophage[n] = len(macrophage)
    QOIs.dendritic[n] = len(dendritic)
    QOIs.CD8[n] = len(CD8)
    QOIs.CD4[n] = len(CD4)

def plotStrain(mcds,n,QOIs,figure,axes,filenameOut,SavePNG):
    Lcell_size = 10;
    Dcell_size = 5;
    Ecell_size = 12;
    RadiusSize = 400.0

    cx = mcds.data['discrete_cells']['position_x'];
    cy = mcds.data['discrete_cells']['position_y'];
    cycle = mcds.data['discrete_cells']['cycle_model']
    current_phase = mcds.data['discrete_cells']['current_phase']
    current_phase = current_phase.astype(int)
    elapsed_time_in_phase = mcds.data['discrete_cells']['elapsed_time_in_phase']
    cell_type = mcds.data['discrete_cells']['cell_type']
    cell_type = cell_type.astype(int)

    strain = mcds.data['discrete_cells']['mechanical_strain']
    pressure = mcds.data['discrete_cells']['simple_pressure']

    melanoma = np.argwhere( (cell_type==2) & (cycle < 100) ).flatten()
    lung = np.argwhere( (cell_type==1) & (cycle < 100) ).flatten()


    QOIs.lung[n] = len(lung)
    QOIs.melanoma[n] = len(melanoma)

    plt.close()
    figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(8,8))
    ax = plt.subplot(111)
    sc = ax.scatter(cx[melanoma],cy[melanoma],c=strain[melanoma], vmin=strain.min(), vmax=strain.max(), s=5, cmap="jet",zorder=1)
    #sc = ax.scatter(cx[lung],cy[lung],c=strain[lung], vmin=strain.min(), vmax=strain.max(), s=5, cmap="jet",zorder=2)
    ax.set_aspect('equal', adjustable='box')
    plt.colorbar(sc, label="strain ($\mu m$)")

    figure, axes = plt.subplots(nrows=1, ncols=1,figsize=(8,8))
    ax = plt.subplot(111)
    sc = ax.scatter(cx[melanoma],cy[melanoma],c=pressure[melanoma], vmin=pressure.min(), vmax=pressure.max(), s=5, cmap="jet",zorder=1)
    #sc = ax.scatter(cx[lung],cy[lung],c=pressure[lung], vmin=pressure.min(), vmax=pressure.max(), s=5, cmap="jet",zorder=2)
    ax.set_aspect('equal', adjustable='box')
    plt.colorbar(sc, label="simple pressure")

    plt.xlim(-RadiusSize, RadiusSize)
    plt.ylim(-RadiusSize, RadiusSize)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    if (SavePNG):
      figure.savefig(filenameOut)
      #plt.savefig(filenameOut)
    else:
      plt.show()
    plt.clf()


def plotAll(mcds,n,QOIs,figure,axes,filenameOut,SavePNG):
    Lcell_size = 10;
    Dcell_size = 5;
    Ecell_size = 12;
    RadiusSize = 400.0

    cx = mcds.data['discrete_cells']['position_x'];
    cy = mcds.data['discrete_cells']['position_y'];
    cycle = mcds.data['discrete_cells']['cycle_model']
    current_phase = mcds.data['discrete_cells']['current_phase']
    current_phase = current_phase.astype(int)
    elapsed_time_in_phase = mcds.data['discrete_cells']['elapsed_time_in_phase']
    cell_type = mcds.data['discrete_cells']['cell_type']
    cell_type = cell_type.astype(int)

    melanoma = np.argwhere( (cell_type==2) & (cycle < 100) ).flatten()
    lung = np.argwhere( (cell_type==1) & (cycle < 100) ).flatten()
    dead = np.argwhere( (cycle >= 100) & ((cell_type==1) | (cell_type==2)) ).flatten()
    macrophage = np.argwhere( cell_type==4 ).flatten()
    dendritic = np.argwhere( cell_type==5 ).flatten()
    CD8 = np.argwhere( cell_type==3 ).flatten()
    CD4 = np.argwhere( cell_type==6 ).flatten()

    QOIs.lung[n] = len(lung)
    QOIs.melanoma[n] = len(melanoma)
    QOIs.dead[n] = len(dead)
    QOIs.macrophage[n] = len(macrophage)
    QOIs.dendritic[n] = len(dendritic)
    QOIs.CD8[n] = len(CD8)
    QOIs.CD4[n] = len(CD4)


    figure.suptitle( '#Lung:'+str("%04i"%(QOIs.lung[n]))+'  #Cancer:'+str("%04i"%(QOIs.melanoma[n]))+ '--  #M:'+str("%04i"%(QOIs.macrophage[n]))+ '  #DC:'+str("%04i"%(QOIs.dendritic[n]))+'\n'+'  #CD8:'+str("%04i"%(QOIs.CD8[n]))+ '  #CD4:'+str("%04i"%(QOIs.CD4[n]))+ '  #Dead:'+str("%04i"%(QOIs.dead[n]))+'  Time:' +str("%8.2f"%(QOIs.times[n]/60.0)) + ' hours', size=14)
    # figure.suptitle( '#NC:'+str("%04i"%(lung_count[n]))+'  #MC:'+str("%04i"%(melanoma_count[n]))+ '--  #M:'+str("%04i"%(macrophage_count[n]))+ '  #DC:'+str("%04i"%(dendritic_count[n]))+'\n'+'  #CD8:'+str("%04i"%(CD8_count[n]))+ '  #CD4:'+str("%04i"%(CD4_count[n]))+ '  #D:'+str("%04i"%(dead_count[n]))+'  Time:' +str("%8.2f"%(n)) + ' hours', size=14)

    ax = plt.subplot(121)
    plt.scatter( cx[lung],cy[lung],c='blue',s=Lcell_size,label='NC');
    plt.scatter( cx[melanoma],cy[melanoma],c='yellow',s=Lcell_size,label='MC');
    plt.scatter( cx[dead],cy[dead],c='black',s=Dcell_size,label='D', alpha=0.25);
    plt.scatter( cx[macrophage],cy[macrophage],c='green',s=Ecell_size,label='M' );
    plt.scatter( cx[dendritic],cy[dendritic],c='brown',s=Ecell_size,label='DC' );
    plt.scatter( cx[CD8],cy[CD8],c='red',s=Ecell_size,label='CD8' );
    plt.scatter( cx[CD4],cy[CD4],c='orange',s=Ecell_size,label='CD4' );
    plt.xlim(-RadiusSize, RadiusSize)
    plt.ylim(-RadiusSize, RadiusSize)
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
    ax.set_aspect('equal')
    plt.subplots_adjust(left=0.08,right=0.93,bottom=0.06,top=0.89,wspace=0.36,hspace=0.26)

    plt.subplot(222)
    TNF = mcds.get_concentrations( 'TNF' );
    X1,Y1 = mcds.get_2D_mesh();

    # v1 = np.linspace(0, TNF.max(), 10, endpoint=True)
    v1 = np.linspace(0, 1.0, 10, endpoint=True)
    plt.contourf(X1,Y1,TNF[:,:,0],levels=v1,cmap='binary');
    # cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.2f"))
    cbar = plt.colorbar(ticks=v1)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.update_ticks()
    plt.xlim(-RadiusSize, RadiusSize)
    plt.ylim(-RadiusSize, RadiusSize)
    plt.title("TNF")

    plt.subplot(224)
    debris = mcds.get_concentrations( 'debris' );
    v1 = np.linspace(0, 1.0, 10, endpoint=True)
    #v1 = np.linspace(0, debris.max(), 10, endpoint=True)
    plt.contourf(X1,Y1,debris[:,:,0],v1,cmap='binary');
    # cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.2f"))
    cbar = plt.colorbar(ticks=v1)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.update_ticks()
    plt.xlim(-RadiusSize, RadiusSize)
    plt.ylim(-RadiusSize, RadiusSize)
    plt.title("Debris")

    if (SavePNG):
      figure.savefig(filenameOut)
      #plt.savefig(filenameOut)
    else:
      plt.draw()
      plt.waitforbuttonpress(0) # this will wait for indefinite time
      #plt.pause(0.2)
    plt.clf()

def plotAll_hist(mcds,n,QOIs,figure,axes,filenameOut,SavePNG):
    Lcell_size = 10;
    Dcell_size = 5;
    Ecell_size = 12;
    RadiusSize = 400.0

    cx = mcds.data['discrete_cells']['position_x'];
    cy = mcds.data['discrete_cells']['position_y'];
    cycle = mcds.data['discrete_cells']['cycle_model']
    current_phase = mcds.data['discrete_cells']['current_phase']
    current_phase = current_phase.astype(int)
    elapsed_time_in_phase = mcds.data['discrete_cells']['elapsed_time_in_phase']
    cell_type = mcds.data['discrete_cells']['cell_type']
    cell_type = cell_type.astype(int)

    melanoma = np.argwhere( (cell_type==2) & (cycle < 100) ).flatten()
    lung = np.argwhere( (cell_type==1) & (cycle < 100) ).flatten()
    dead = np.argwhere( (cycle >= 100) & ((cell_type==1) | (cell_type==2)) ).flatten()
    macrophage = np.argwhere( cell_type==4 ).flatten()
    dendritic = np.argwhere( cell_type==6 ).flatten()
    CD8 = np.argwhere( cell_type==3 ).flatten()
    CD4 = np.argwhere( cell_type==7 ).flatten()

    QOIs.lung[n] = len(lung)
    QOIs.melanoma[n] = len(melanoma)
    QOIs.dead[n] = len(dead)
    QOIs.macrophage[n] = len(macrophage)
    QOIs.dendritic[n] = len(dendritic)
    QOIs.CD8[n] = len(CD8)
    QOIs.CD4[n] = len(CD4)


    figure.suptitle( '#NC:'+str("%04i"%(QOIs.lung[n]))+'  #MC:'+str("%04i"%(QOIs.melanoma[n]))+ '--  #M:'+str("%04i"%(QOIs.macrophage[n]))+ '  #DC:'+str("%04i"%(QOIs.dendritic[n]))+'\n'+'  #CD8:'+str("%04i"%(QOIs.CD8[n]))+ '  #D:'+str("%04i"%(QOIs.dead[n]))+'  Time:' +str("%8.2f"%(QOIs.times[n]/60.0)) + ' hours', size=14)
    # figure.suptitle( '#NC:'+str("%04i"%(lung_count[n]))+'  #MC:'+str("%04i"%(melanoma_count[n]))+ '--  #M:'+str("%04i"%(macrophage_count[n]))+ '  #DC:'+str("%04i"%(dendritic_count[n]))+'\n'+'  #CD8:'+str("%04i"%(CD8_count[n]))+ '  #CD4:'+str("%04i"%(CD4_count[n]))+ '  #D:'+str("%04i"%(dead_count[n]))+'  Time:' +str("%8.2f"%(n)) + ' hours', size=14)

    ax = plt.subplot2grid((4,4), (0,0), colspan=3, rowspan=3)
    plt.scatter( cx[lung],cy[lung],c='blue',s=Lcell_size,label='NC');
    plt.scatter( cx[melanoma],cy[melanoma],c='yellow',s=Lcell_size,label='MC');
    plt.scatter( cx[dead],cy[dead],c='black',s=Dcell_size,label='D', alpha=0.25);
    plt.scatter( cx[macrophage],cy[macrophage],c='green',s=Ecell_size,label='M' );
    plt.scatter( cx[dendritic],cy[dendritic],c='brown',s=Ecell_size,label='DC' );
    plt.scatter( cx[CD8],cy[CD8],c='red',s=Ecell_size,label='CD8' );
    # plt.scatter( cx[CD4],cy[CD4],c='orange',s=Ecell_size,label='CD4' );
    plt.xlim(-RadiusSize, RadiusSize)
    plt.ylim(-RadiusSize, RadiusSize)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.legend(bbox_to_anchor =(0.75, 1.15), ncol = 4)
    ax.set_aspect('equal')
    plt.subplots_adjust(left=0.08,right=0.93,bottom=0.06,top=0.84,wspace=0.36,hspace=0.45)

    ax=plt.subplot2grid((4,4), (0,3))
    Neo = mcds.get_concentrations( 'neoantigens' );
    X1,Y1 = mcds.get_2D_mesh();
    if (Neo.max() > 0):
      v1 = np.linspace(0, Neo.max(), 10, endpoint=True)
      # v1 = np.linspace(0, 1.0, 10, endpoint=True)
      plt.contourf(X1,Y1,Neo[:,:,0],levels=v1,cmap='binary');
      # cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.2f"))
      cbar = plt.colorbar(ticks=v1)
      cbar.formatter.set_powerlimits((0, 0))
      cbar.update_ticks()
    else:
      plt.contourf(X1,Y1,Neo[:,:,0],cmap='binary');
      # cbar = plt.colorbar(format=FormatScalarFormatter("%.2f"))
      cbar = plt.colorbar()
    plt.xlim(-RadiusSize, RadiusSize)
    plt.ylim(-RadiusSize, RadiusSize)
    plt.title("Neoantigens")
    ax.set_aspect('equal')

    ax=plt.subplot2grid((4,4), (1,3))
    chemokine = mcds.get_concentrations( 'chemokine' );
    X1,Y1 = mcds.get_2D_mesh();
    if (chemokine.max() > 0):
      v1 = np.linspace(0, chemokine.max(), 10, endpoint=True)
      # v1 = np.linspace(0, 1.0, 10, endpoint=True)
      plt.contourf(X1,Y1,chemokine[:,:,0],levels=v1,cmap='binary');
      # cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.2f"))
      cbar = plt.colorbar(ticks=v1)
      cbar.formatter.set_powerlimits((0, 0))
      cbar.update_ticks()
    else:
      plt.contourf(X1,Y1,chemokine[:,:,0],cmap='binary');
      # plt.colorbar(format=FormatScalarFormatter("%.2f"))
      cbar = plt.colorbar()
    plt.xlim(-RadiusSize, RadiusSize)
    plt.ylim(-RadiusSize, RadiusSize)
    plt.title("Chemokine")
    ax.set_aspect('equal')

    ax=plt.subplot2grid((4,4), (2,3))
    PIC = mcds.get_concentrations( 'pro-inflammatory cytokine' );
    if (PIC.max() > 0):
      v1 = np.linspace(0, PIC.max(), 10, endpoint=True)
      # v1 = np.linspace(0, 1.0, 10, endpoint=True)
      plt.contourf(X1,Y1,PIC[:,:,0],v1,cmap='binary');
      # cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.2f"))
      cbar = plt.colorbar(ticks=v1)
      cbar.formatter.set_powerlimits((0, 0))
      cbar.update_ticks()
    else:
      plt.contourf(X1,Y1,PIC[:,:,0],cmap='binary');
      # cbar = plt.colorbar(format=FormatScalarFormatter("%.2f"))
      cbar = plt.colorbar()
    plt.xlim(-RadiusSize, RadiusSize)
    plt.ylim(-RadiusSize, RadiusSize)
    plt.title("Pro-inf cytokine")
    ax.set_aspect('equal')

    ax=plt.subplot2grid((4,4), (3,0))
    debris = mcds.get_concentrations( 'debris' );
    if (debris.max() > 0):
      #v1 = np.linspace(0, 1.0, 10, endpoint=True)
      v1 = np.linspace(0, debris.max(), 10, endpoint=True)
      plt.contourf(X1,Y1,debris[:,:,0],v1,cmap='binary');
      # cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.2f"))
      cbar = plt.colorbar(ticks=v1)
      cbar.formatter.set_powerlimits((0, 0))
      cbar.update_ticks()
    else:
      plt.contourf(X1,Y1,debris[:,:,0],cmap='binary');
      # cbar = plt.colorbar(format=FormatScalarFormatter("%.2f"))
      cbar = plt.colorbar()
    plt.xlim(-RadiusSize, RadiusSize)
    plt.ylim(-RadiusSize, RadiusSize)
    plt.title("Debris")
    ax.set_aspect('equal')

    ax=plt.subplot2grid((4,4), (3,1))
    AIC = mcds.get_concentrations( 'anti-inflammatory cytokine' );
    if (AIC.max() > 0):
      v1 = np.linspace(0, AIC.max(), 10, endpoint=True)
      # v1 = np.linspace(0, 1.0, 10, endpoint=True)
      plt.contourf(X1,Y1,AIC[:,:,0],v1,cmap='binary');
      # cbar = plt.colorbar(ticks=v1,format=FormatScalarFormatter("%.2f"))
      cbar = plt.colorbar(ticks=v1)
      cbar.formatter.set_powerlimits((0, 0))
      cbar.update_ticks()
    else:
      plt.contourf(X1,Y1,AIC[:,:,0],cmap='binary');
      # cbar = plt.colorbar(format=FormatScalarFormatter("%.2f"))
      cbar = plt.colorbar()
      cbar.formatter.set_powerlimits((0, 0))
      cbar.update_ticks()
    plt.xlim(-RadiusSize, RadiusSize)
    plt.ylim(-RadiusSize, RadiusSize)
    plt.title("Anti-inf cytokine")
    ax.set_aspect('equal')

    plt.subplot2grid((4,4), (3,2))
    sns.set()
    sns.set_style('white')
    plt.xlim(0, 1)
    sns.distplot(NeoCell[melanoma],color='black').get_lines()[0]
    plt.title("Neoantigen")

    plt.subplot2grid((4,4), (3,3))
    plt.xlim(0, 1)
    if (np.std(PDL1[melanoma]) < 1.0e-15):
        print("\n")
    else:
        print(np.std(PDL1[melanoma]))
        sns.distplot(PDL1[melanoma],color='black').get_lines()[0]
    plt.title("PDL1 expression")

    if (SavePNG):
      figure.savefig(filenameOut)
      #plt.savefig(filenameOut)
    else:
      plt.draw()
      plt.waitforbuttonpress(0) # this will wait for indefinite time
      #plt.pause(0.2)
    plt.clf()

def plotCurves(QOIs,STD):
    #fig, ax = plt.subplots()
    # plt.plot(QOIs.times/1440.0,QOIs.lung, color='blue',label='normal')
    # plt.fill_between(QOIs.times/1440.0, QOIs.lung - STD.lung, QOIs.lung + STD.lung, alpha=0.2, color='blue')
    plt.plot(QOIs.times/1440.0,QOIs.melanoma, color='gold',label='melanoma')
    plt.fill_between(QOIs.times/1440.0, QOIs.melanoma - STD.melanoma, QOIs.melanoma + STD.melanoma, alpha=0.2, color='gold')
    plt.plot(QOIs.times/1440.0,QOIs.dead, color='black',label='dead')
    plt.fill_between(QOIs.times/1440.0, QOIs.dead - STD.dead, QOIs.dead + STD.dead, alpha=0.2, color='black')
    plt.plot(QOIs.times/1440.0,QOIs.macrophage, color='green',label='macrophage')
    plt.fill_between(QOIs.times/1440.0, QOIs.macrophage - STD.macrophage, QOIs.macrophage + STD.macrophage, alpha=0.2, color='green')
    plt.plot(QOIs.times/1440.0,QOIs.dendritic, color='brown',label='dendritic')
    plt.fill_between(QOIs.times/1440.0, QOIs.dendritic - STD.dendritic, QOIs.dendritic + STD.dendritic, alpha=0.2, color='brown')
    plt.plot(QOIs.times/1440.0,QOIs.CD8, color='red',label='CD8')
    plt.fill_between(QOIs.times/1440.0, QOIs.CD8 - STD.CD8, QOIs.CD8 + STD.CD8, alpha=0.2, color='red')
    plt.plot(QOIs.times/1440.0,QOIs.CD4, color='orange',label='CD4')
    plt.fill_between(QOIs.times/1440.0, QOIs.CD4 - STD.CD4, QOIs.CD4 + STD.CD4, alpha=0.2, color='orange')
    plt.xlabel("Time (days)")
    plt.ylabel("Number of cell")
    plt.legend(loc='upper left')
    plt.show()

def plotExternalImmune():
    filename= folder+'/dm_tc.dat'
    input = np.loadtxt(filename, dtype='f', delimiter=' ')
    col1 = np.array(input[:,0])
    col2 = np.array(input[:,1])
    col3 = np.array(input[:,1])
    col4 = np.array(input[:,1])
    col5 = np.array(input[:,1])
    col6 = np.array(input[:,1])
    time = np.arange(len(col1))
    plt.plot(time,col1,label='DM')
    plt.plot(time,col2,label='TC')
    plt.plot(time,col3,label='TH1')
    plt.plot(time,col4,label='TH2')
    plt.plot(time,col5,label='TCt')
    plt.plot(time,col6,label='Tht')
    plt.xlabel("Time (hours)")
    plt.ylabel("Number of cell")
    plt.legend(loc='upper left')
    plt.show()

def Read_files(initial_index,last_index,folder,SavePNG,func=plotAll):
    QOIs = NumberOfCells(last_index-initial_index+1)
    figure, axes = plt.subplots(nrows=2, ncols=2,figsize=(10,8))
    for n in range( initial_index,last_index+1 ):
        filenameOut=folder+'/output'+"%08i"%n+'.png'
        filename= "output%08i"%n+'.xml'
        mcds=pyMCDS(filename,folder)
        QOIs.times[n]= mcds.get_time()
        func(mcds,n,QOIs,figure,axes,filenameOut,SavePNG)

if __name__ == '__main__':
    if (len(sys.argv) != 5):
      print("Please provide 4 args: InitialTime LastTime, SavePNG (bool), and folder")
      sys.exit(1)
    initial_index = int(sys.argv[1]);
    last_index = int(sys.argv[2]);
    SavePNG = int(sys.argv[3])
    folder = sys.argv[4]

    Read_files(initial_index,last_index,folder,SavePNG)
    #Read_files(initial_index,last_index,folder,SavePNG,func=plotAll_hist)
    #Read_files(initial_index,last_index,folder,SavePNG,func=plotStrain)
    #Neoantigens_plot(initial_index,last_index,folder)
    #Replicates_plot(initial_index,last_index,folder,5)
