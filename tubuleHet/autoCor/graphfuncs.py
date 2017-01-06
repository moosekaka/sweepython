# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 01:22:22 2015

@author: sweel
"""

def getBptsEpts(vtkData, curGrph):
    """Return branchpoints and endpoints Indexes

    Parameters
    ----------
    vtkData :
        Vtk polydata cell
    curGrph :
        Network x graph
    """
    bpts = [
        curGrph.node[i]['coord'] for i in curGrph.nodes()
        if curGrph.node[i]['degree'] > 1]
    bptsId = [vtkData.find_point(el) for el in bpts]

    epts = [
        curGrph.node[i]['coord'] for i in curGrph.nodes()
        if curGrph.node[i]['degree'] == 1]
    eptsId = [vtkData.find_point(el) for el in epts]

    return(bptsId, eptsId)
