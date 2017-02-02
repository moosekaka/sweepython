# -*- coding: utf-8 -*-
"""
Module for functions for mito network normalization
@author: sweel_rafelski
"""
import os
import os.path as op
import vtk
import vtk.util.numpy_support as vnpy
from numpy import ceil, mean, percentile
# pylint: disable=C0103
datadir = op.join(os.getcwd())


def vtk_read(fpath, readertype='vtkPolyDataReader'):
    """
    Reads vtk files and returns the VTK object
    """
    reader = getattr(vtk, readertype)()
    reader.SetFileName(fpath)
    reader.Update()
    data = reader.GetOutput()
    return data


def point_cloud_scalars(skelpath, ch1path, ch2path, **kwargs):
    """
    Returns scalar values from voxels data (eg. *resampledVTK*) lying within
    a point cloud of a specified radius for each point

    Parameters
    ----------
    skelpath, ch1path, ch2path : str
        filepaths to skeleton and volume/voxels VTK output of respective
        channels to be normalized
    kwargs :
        radius value argument (float) for _pointcloud()

        *will default use a value of 2.5 pixels*
    Returns
    -------
    polydata : VTK poly
        polydata object with raw scalar values and width
    np_voxel1, np_voxel2 : Numpy Array
        numpy array conversion of VTK voxel intensities
    """
    dataSkel = vtk_read(skelpath)
    voxels_ch1 = vtk_read(ch1path,
                          readertype='vtkStructuredPointsReader')
    voxels_ch2 = vtk_read(ch2path,
                          readertype='vtkStructuredPointsReader')

    ptcld_ch1, ptcld_ch2 = _pointcloud(dataSkel, voxels_ch1, voxels_ch2,
                                       **kwargs)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(dataSkel.GetPoints())
    polydata.SetLines(dataSkel.GetLines())
    polydata.GetPointData().AddArray(ptcld_ch1)
    polydata.GetPointData().AddArray(ptcld_ch2)
    np_voxel1 = vnpy.vtk_to_numpy(voxels_ch1.GetPointData().GetScalars())
    np_voxel2 = vnpy.vtk_to_numpy(voxels_ch2.GetPointData().GetScalars())
    polydata.GetPointData().AddArray(dataSkel.
                                     GetPointData().GetScalars('Width'))
    polydata.GetPointData().GetArray(2).SetName("tubewidth")
    return polydata, np_voxel1, np_voxel2


def _pointcloud(skel, ch1, ch2, radius=2.5):
    vox_ch2 = vtk.vtkDoubleArray()
    vox_ch2.SetName("vox_ch2")
    vox_ch1 = vtk.vtkDoubleArray()
    vox_ch1.SetName("vox_ch1")

    # vtk.pointlocator() used to to find the set of points lying within the
    # radius of the point of interest, returns results in a list called
    # result
    loc = vtk.vtkPointLocator()
    loc.SetDataSet(ch2)
    loc.BuildLocator()

    inten_ch1 = ch1.GetPointData().GetScalars().GetTuple1
    inten_ch2 = ch2.GetPointData().GetScalars().GetTuple1
    result = vtk.vtkIdList()

#   averaging of pts intensity value surrounding each point in skel
    for n in range(skel.GetNumberOfPoints()):

        pt_of_int = tuple(ceil(i/.055) for i in skel.GetPoint(n))
        loc.FindPointsWithinRadius(radius, pt_of_int, result)
        vox_id = [result.GetId(i) for i in range(result.GetNumberOfIds())]

        vox_ch1.InsertNextValue(mean([inten_ch1(m) for m in vox_id]))

        vox_ch2.InsertNextValue(mean([inten_ch2(m) for m in vox_id]))
    return vox_ch1, vox_ch2


def normalize_skel(polydata, raw_vox_ch1, raw_vox_ch2,
                   background_thresh=5., **kwargs):
    """
    Normalize channels to correct for focal plane intensity variations

    Parameters
    ----------
    polydata : vtkPolyData
        vtk object returned from pt_cld_sclrs()
    raw_vox_ch1, raw_vox_ch2 : Numpy array
        Voxel intensity values in numpy array format
    backgrnd_thresh : float
        The default threshold of a background value is the 5th percentile of
        voxel intensities in the respective channel. This might have to be
        changed depending on the experimental conditions

    kwargs
    ------
    background: Str
        path to background file provided by user, should be a dictionary with
        channel labels as keys
    """
    temp = polydata.GetPointData()
    vox_ch1 = vnpy.vtk_to_numpy(temp.GetArray('vox_ch1'))
    vox_ch2 = vnpy.vtk_to_numpy(temp.GetArray('vox_ch2'))

    # Check to see if background file provided by user
    background = kwargs.pop('backgroundfile')
    if background:
        if background['ch2'] > min(vox_ch2):
            min_ch2 = min(vox_ch2)-1.  # ensure no zero divisions
        else:
            min_ch2 = background['ch2']-1.  # USUAL CASE

        min_ch1 = min(background['ch1'], min(vox_ch1))
        print ("Background value of gfp={:4.0f} & rfp={:4.0f} was used".
               format(background['ch1'], background['ch2']))

    else:
        min_ch1 = percentile(raw_vox_ch1, background_thresh)
        min_ch2 = percentile(raw_vox_ch2, background_thresh)

    # background Substracted rfp and gfps
    ch2_bckgrnd = vox_ch2 - min_ch2
    ch1_bckgrnd = vox_ch1 - min_ch1

    # width equivalent
    width_eqv = ch2_bckgrnd / min(ch2_bckgrnd)
    unscaled_dy = ch1_bckgrnd / width_eqv  # raw DY/W normalized values

    # rescale DY to minmax
    _min = min(unscaled_dy)
    _max = max(unscaled_dy)
    normalized_dy = ((unscaled_dy - _min)/(_max - _min))
    tubewidth = temp.GetArray('tubewidth')

    # Output results as a labelled dictionary
    results = {'normalized_dy': normalized_dy,
               'unscaled_dy': unscaled_dy,
               'ch1_bckgrnd': ch1_bckgrnd,
               'ch2_bckgrnd': ch2_bckgrnd,
               'width_eqv': width_eqv,
               'tubewidth': tubewidth}
    return results


def write_vtk(dat, fname, **kwargs):
    """
    Write out a vtk file using VTK polydata object *dat* and a filename *fname*
    with optional labels dictionary *kwargs* for the outputs

    kwargs
    ------
    Default dictionary keys are:

    * normalized_dy
    * unscaled_dy
    * ch1_bckgrnd
    * ch2_bckgrnd
    * width_eqv'

    """
    for k in sorted(kwargs):
        try:
            temp = vnpy.numpy_to_vtk(kwargs[k])
        except TypeError:  # already a VTK type
            temp = kwargs[k]
        temp.SetName(k)
        dat.GetPointData().AddArray(temp)

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(fname)
    writer.SetInputData(dat)
    writer.Update()
