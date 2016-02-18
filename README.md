# What is in this repository?
**A set of open-source computational tools developed in Python to perform *multi-scale* analysis of *structure-function* relationship in mitochondrial networks.**

##Segmentation of mitochondrial networks from live yeast cells 3D images taken with a *spinning disk confocal* microscope.
- This segmentation utilizes [`MitoGraph`](https://github.com/vianamp/MitoGraph.git), a C++ skeletonization and segmentation program developed by Matheus Viana. It has been fully validated in yeast cells \([Viana, Lim *et al.*](http://www.ncbi.nlm.nih.gov/pubmed/25640425)). 

- The raw data from MitoGraph was 'munged' into a database using [`Pandas`](http://pandas.pydata.org/pandas-docs/stable/) before further calculations of various statistical and functional parameters.

###Modulation of mitochondrial function by altering carbon source growth conditions.
- We grew cells in growth conditions utilizing different carbon sources in order to modulate their respiratory conditions [`O2 results`](https://github.com/moosekaka/sweepython/tree/master/o2).

###Investigating the link between structure of mitochondria and heterogeneity of mitochondrial membrane potential (Δψ).
- We showed the non-random heterogenous distribution of Δψ within a single mitochondrial tubule [`example here`](https://github.com/moosekaka/sweepython/tree/master/tubuleHet).

- We also detail our investigation into the relationship beteen heterogeneity of mitochondrial function and network topology [`example here`](https://github.com/moosekaka/sweepython/tree/master/networkHet).

###Analysis of mitochondrial functional asymmetry in mother and daughter yeast cells.
- A [`set of tools`](https://github.com/moosekaka/sweepython/tree/master/vtk_viz) to interactively visualize the 3D mitochondrial skeleton and pick points to classify the cell as a mother or daughter region. This [`Ipython notebook`](https://github.com/moosekaka/sweepython/blob/master/mother%20bud%20analysis.ipynb) demonstrates an application of these tools to study how mitochondrial membrane potential (Δψ) is distributed differently between the mother and daughter cell.  

<p align="center">
<img src="https://github.com/moosekaka/sweepython/blob/master/images/coors.png" width="300" />
</p>

####Requirements
* Pandas
* Matplotlib
* Mayavi\**
* NetworkX
* Numpy
* Seaborn
* Scipy
* VTK

_\** MayaVi has a known conflict when run under IPython and Python 2.7, specifically incompatible API versions. To fix, set QT_API=pyqt under your enviroment variables, either in BASH or Windows advanced settings._


The best way to ensure all dependencies are fulfilled is by installing the [Anaconda](https://www.continuum.io/downloads) Python package.

