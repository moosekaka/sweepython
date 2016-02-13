# What is in this repository?
**A set of open-source computational tools developed in Python to perform *multi-scale* analysis of *structure-function* relationship in mitochondrial networks.**

##Segmentation of mitochondrial networks from live yeast cells 3D images taken with a *spinning disk confocal* microscope.
This segmentation utilizes [MitoGraph](https://github.com/vianamp/MitoGraph.git), a C++ skeletonization and segmentation program developed by Matheus Viana. It has been fully validated in yeast cells \([*Viana, Lim et al*](http://www.ncbi.nlm.nih.gov/pubmed/25640425)\). 

##Data Munging in Pandas for data exploration and visualization

###Modulation of mitochondrial function by altering carbon source growth conditions.

###Investigating tubule scale heterogeneity of mitochondrial membrane potential (Δψ).

###Analysis of mitochondrial functional asymmetry in mother and daughter yeast cells.
A [set of tools](https://github.com/moosekaka/sweepython/tree/master/vtk_viz) to interactively visualize the 3D mitochondrial skeleton and pick points to classify the cell as a mother or daughter region.

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

