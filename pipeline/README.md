## A set of tools used for a pipeline to map function (mito. membrane potential ΔΨ) to structure
- This notebook shows a [`Python pipeline`](https://github.com/moosekaka/sweepython/blob/master/Pipeline.ipynb) to map functional values onto the segmented skeleton, with corrections for optical defects.
- This pipeline utilizes Mayavi's wrapped version of VTK to take advantage of Pythonic coding syntax when working with VTK data.

### Preprocessing
- The data are 3d stacks of yeast cell images taken with a spinning disk confocal microscope. The cells express a fluorescent protein that is targeted to the matrix of the mitochondria. They are also labelled with a ΔΨ potential sensitive dye.
- Individual cell ROI's were picked using [`FIJI`](http://fiji.sc/Fiji)
- The input files ('\*.vtk') were generated using [`MitoGraph v2.0`](http://www.rafelski.com/susanne/MitoGraph.html).

### Pipeline details
- The core functions are in `modules pipefuncs.py` and `make_network.py`. These function perform point cloud averaging to assign scalar values to the skeleton, as well as optical correction and normalization to get an accurate 'heatmap' of the distribution of ΔΨ along the skeleton.
- In addition for each skeleton a graph network is constructed by examining every endpoints of every linesegment in the skeleton for coincident points. A vectorized implementation using [`Numpy's broadcasting`](http://scipy.github.io/old-wiki/pages/EricsBroadcastingDoc) feature speeds up this implementation for large networks.
