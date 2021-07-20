# A multiscale model of melanoma micrometastases with immune system interaction
Building off of recent COVID-19 modeling (https://www.biorxiv.org/content/10.1101/2020.04.02.019075v4), we developed an agent-based model of melanoma metastases in lung tissue and immune response. The model includes cancer cells, lung cells, key immune cell types, and cytokine-mediated and mechanical interactions. Immune cells traffic to and from the lymphatic system to drive an expanding immune response. 

This model is part of the pilot project for building a Cancer patient digital twins (CPDTs) framework.

![alt ensure executable](https://raw.githubusercontent.com/heberlr/melanoma/master/beta/model_scheme.png)

To compile and run [PhysiCell](http://physicell.mathcancer.org/) for this model:
```
# your compiler needs to support OpenMP
$ make
$ ./melanoma
```