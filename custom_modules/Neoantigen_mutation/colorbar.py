import numpy as np
import matplotlib.pyplot as plt

Data = []
for line in open('Data256.dat', 'r'):
  values = [float(s) for s in line.split()]
  Data.append(values)
Data = np.array(Data)

Size = np.unique(Data[:,0])
Xtick = np.unique(Data[:,2])
Ytick = np.unique(Data[:,1])
DataMatrix = np.reshape(Data[:,3], (Xtick.size,Ytick.size))

fig, ax = plt.subplots()
im = ax.imshow(DataMatrix,cmap='winter',origin = 'lower')

# We want to show all ticks...
ax.set_xticks(np.arange(Xtick.size))
ax.set_yticks(np.arange(Ytick.size))
# ... and label them with the respective list entries
ax.set_xticklabels(Xtick)
ax.set_yticklabels(Ytick)
#
# # Rotate the tick labels and set their alignment.
# plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#          rotation_mode="anchor")
#
#Loop over data dimensions and create text annotations.
for i in range(Xtick.size):
    for j in range(Ytick.size):
        text = ax.text(j, i, DataMatrix[i, j],
                       ha="center", va="center", color="white")

v1 = np.linspace(0, DataMatrix.max(), 9, endpoint=True)
plt.colorbar(im, ticks=v1, label='Hamming distance')
plt.title( 'Case with '+str("%03i"%(Size))+' elements' )
plt.xlabel('Generation')
plt.ylabel('Mutation rate')
fig.tight_layout()
plt.show()
