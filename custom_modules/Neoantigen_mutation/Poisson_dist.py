import numpy as np
from scipy.special import factorial
import matplotlib.pyplot as plt
import seaborn as sns

Data = []
for line in open('PoissonDist.dat', 'r'):
  values = [float(s) for s in line.split()]
  Data.append(values)
Data = np.hstack(Data)
prob_mutation = 10.0

print(Data.shape)

fig, ax = plt.subplots(1,2, figsize=(10,4))

sns.distplot(Data,color='green',ax=ax[0], kde=False, hist=True)
x = np.arange(16)
PoissonResut = np.exp(-0.5*prob_mutation)*np.power(0.5*prob_mutation,x)/factorial(x)
print(PoissonResut)
ax[1].plot(x,PoissonResut)
plt.show()
