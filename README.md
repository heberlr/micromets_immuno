# A multiscale model of pulmonary micrometastases with immune system interaction

## The multiscale model
Building off of recent COVID-19 modeling (https://www.biorxiv.org/content/10.1101/2020.04.02.019075v4), we developed an agent-based model of micrometastases in lung tissue and immune response. The model includes cancer cells, lung cells, key immune cell types, and cytokine-mediated and mechanical interactions. Immune cells traffic to and from the lymphatic system to drive an expanding immune response.

This model is part of the pilot project for building a Cancer patient digital twins (CPDTs) framework.

![alt ensure executable](https://raw.githubusercontent.com/heberlr/micromets_lung/development/beta/model_scheme.png)

To compile and run [PhysiCell](http://physicell.mathcancer.org/) for this model:

```
# your compiler needs to support OpenMP
$ make
$ ./micromets_lung
```
## Data analysis
**Requirements**: our analysis have been performed using *Python 3.9.18*. The following modules are required to run the [jupyter notebook](https://github.com/heberlr/micromets_lung/blob/development/Data_Analysis/PlotResults.ipynb):
- Matplotlib 3.8.0 (https://matplotlib.org/)
- SciPy 1.11.2 (https://scipy.org/)
- pandas 1.5.3 (https://pandas.pydata.org/)
- wget 3.2 (https://pypi.org/project/wget/)
- torchvision 0.15.2 (https://pytorch.org/vision/stable/index.html)
- fitter 1.5.2 (https://fitter.readthedocs.io/en/latest/)
- matplotlib_venn 0.11.9 (https://pypi.org/project/matplotlib-venn/)
- seaborn 0.12.2 (https://seaborn.pydata.org/)

If you are using Anaconda, you can create the conda environment manually:
```
conda create -n micromets python=3.9
conda activate micromets
pip install jupyter pandas==1.5.3 matplotlib==3.8.0 wget==3.2 torchvision==0.15.2 scipy==1.11.2 fitter==1.5.2 matplotlib_venn==0.11.9 seaborn==0.12.2
```
or you can create from the [micromets.yml](https://github.com/heberlr/micromets_lung/blob/development/Data_Analysis/micromets.yml) file:
```
conda env create -f micromets.yml
```
