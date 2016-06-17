# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:17:17 2016

@author: sweel_Rafelski
"""

     #  DataFrame for mom bin points + first point on bud
            fp = cell[cell['type'] == 'bud'][:1].reset_index(drop=True)
            Xcell = cell.groupby('whole_cell_binpos').DY.mean()
            Xcell[Xcell > cell.neckpos_cellaxis.max()] = np.nan
            Xcell = Xcell.reset_index()
            dfXcell = vf.xcell(Xcell, fp)
            dfXcell['type'] = celltype
            dfmom_fp = dfmom_fp.append(dfXcell)  # DataFrame output