# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 17:54:10 2015
       module to pick mom, bud ,neck positions
"""
import os
import os.path as op
import numpy as np
from mayavi import mlab
import pandas as pd
from tvtk.api import tvtk
from collections import defaultdict
from mombud.functions import vtkvizfuncs as vz
import wrappers as wr
# pylint: disable=C0103
vtkF = defaultdict(dict)
datadir = op.join(os.getcwd(), 'mutants')
rawdir = op.join(os.getcwd(), 'mutants')


def picker_callback(picker_obj):
    """callback function
    """
    global tip, base
    picked = picker_obj.actors

    if momshell.actor.actor._vtk_obj in [o._vtk_obj for o in picked]:
        x, y, z = picker_obj.pick_position
#        basetext.set(text='%4.2f %4.2f %4.2f' % (x, y, z),
#                     position=[x, y, z])
        base = [x, y, z]
        cursorbase.actor.actor.set(position=[x, y, z])

    if budshell.actor.actor._vtk_obj in [o._vtk_obj for o in picked]:
        x, y, z = picker_obj.pick_position
#        tiptext.set(text='%4.2f %4.2f %4.2f' % (x, y, z),
#                    position=[x, y, z])
        tip = [x, y, z]
        cursortip.actor.actor.set(position=[x, y, z])


def picker_callback2(picker_obj):
    """callback function
    """
    global neck
    x, y, z = picker_obj.pick_position
#    necktext.set(text='%4.2f %4.2f %4.2f' % (x, y, z),
#                 position=[x, y, z])
    neck = [x, y, z]
    cursorneck.actor.actor.set(position=[x, y, z])


def picker_callback4(picker_obj):
    """callback function to change z position
    """
    global zpt
    _, _, znew = picker_obj.pick_position
    zpt = [_, _, znew]
    xm, ym, zm = momshell.actor.actor.position
    xd, yd, zd = budshell.actor.actor.position
    momshell.actor.actor.set(position=[xm, ym, znew])
    budshell.actor.actor.set(position=[xd, yd, znew])


def change_interaction():
    """
    Function to change the interaction style and adds observers
    """
    istyle = tvtk.InteractorStyleTerrain()
    iactor = figone.scene.interactor
    iactor.interactor_style = istyle
    istyle.add_observer('KeyPressEvent', callback)


def callback(obj, event):
    """
    Callback function, handels the interaction events
    """

    if event == 'KeyPressEvent':
        # Get the pressed key
        key = obj.GetInteractor().GetKeyCode()
        xn, yn, zn = neck
        xb, yb, zb = base
        xt, yt, zt = tip

        try:
            zpos = zpt[2]
        except NameError:
            zpos = zp

        if key == 'd':
            #  remove old arrow
            for i in figone.children:
                if i.name == "VTK Data (PolyData)":
                    if i.children[0].children[0].name == 'mombudaxis':
                        i.remove()
        #  get orientation of the new arrow
            tr, rot, scale = vz.arrowvect(base, tip, neck)
            arrsource = tvtk.ArrowSource(shaft_radius=.01,
                                         shaft_resolution=18,
                                         tip_length=.15,
                                         tip_radius=.05,
                                         tip_resolution=18)
            transformPD = tvtk.TransformPolyDataFilter()
            transformPD = tvtk.TransformPolyDataFilter(input=arrsource.output,
                                                       transform=tr)
            arrdata = transformPD.output
            mlab.pipeline.surface(arrdata,
                                  figure=figone,
                                  opacity=.7,
                                  name='mombudaxis')

        if key == 'x':
            output = op.join(datadir, '%s.csv' % filekey)
            f = open(output, 'w')
            f.write('%s\n' % filekey)
            f.write('neck,%6.4f,%6.4f,%6.4f\n' % (xn, yn, zn))
            f.write('base,%6.4f,%6.4f,%6.4f\n' % (xb, yb, zb))
            f.write('tip,%6.4f,%6.4f,%6.4f\n' % (xt, yt, zt))
            f.write('centerpt,%6.4f\n' % (zpos))
            f.close()
            print 'results recorded!'

# filelist and graph list
vtkF = wr.ddwalk(op.join(rawdir, 'normalizedVTK'), '*skeleton.vtk',
                 start=5, stop=-13)

filekeys = {item: vtkF[media][item] for media
            in sorted(vtkF.keys()) for item
            in sorted(vtkF[media].keys())}

DataSize = pd.read_table(op.join(datadir, 'Results.txt'))
df = DataSize.ix[:, 1:]
df['cell'] = df.ix[:, 'Label'].apply(lambda x: x.partition(':')[2])
df['vol'] = 4/3 * np.pi * (df.Major*.055/2) * (df.Minor*.055/2) ** 2

# Draw cell using cellplot and edgeplot
if __name__ == "__main__":
    mlab.close(all=True)
    filekey = 'WT_032716_009_RFPstack_061'
    df2 = vz.getelipspar(filekey, df)
    df2 = df2.sort_values(by='vol')
    df2.reset_index(drop=True, inplace=True)
    df2.index = ['bud', 'mom']
    df2['center'] = zip((df2.X - 25)*.055, (225 - df2.Y)*.055)

    figone = mlab.figure(figure=filekey,
                         size=(800, 600),
                         bgcolor=(0.1, 0.1, 0.1))
    vtkobj, tubeout = vz.cellplot(figone, filekeys[filekey], rad=.07)
    figone.scene.disable_render = True

    xmin, xmax, ymin, ymax, zmin, zmax = vtkobj.outputs[0].bounds
    zp = (zmax-zmin)/2
    vz.adjustlut(tubeout)
    momshell, moms = vz.drawelips('mom', df2, zpos=zp)
    budshell, buds = vz.drawelips('bud', df2, zpos=zp)
    moms.set(name='shell')
    buds.set(name='shell')

    basetext = mlab.text3d(0, 0, 0,
                           '%4.2f %4.2f %4.2f' % (0, 0, 0),
                           color=(.0, .25, .8),
                           scale=.2)
    tiptext = mlab.text3d(1, 0, 0,
                          '%4.2f %4.2f %4.2f' % (0, 0, 0),
                          color=(.2, .8, .2),
                          scale=.2)
    necktext = mlab.text3d(1, 1, 0,
                           '%4.2f %4.2f %4.2f' % (0, 0, 0),
                           color=(.9, .1, .1),
                           scale=.2)
    cursorneck = mlab.points3d(0, 0, 0,
                               mode='2dcross',
                               scale_factor=.25,
                               color=(.9, .1, .1))
    cursorbase = mlab.points3d(0, 0, 0,
                               mode='2dcross',
                               scale_factor=.25,
                               color=(.0, .25, .9))
    cursortip = mlab.points3d(0, 0, 0,
                              mode='2dcross',
                              scale_factor=.25,
                              color=(.2, .7, .2))

    figone.scene.disable_render = False
    mlab.view(0, 0, 180)
    mlab.view(distance='auto')
    change_interaction()
    figone.on_mouse_pick(picker_callback)
    figone.on_mouse_pick(picker_callback2, button='Right')
    figone.on_mouse_pick(picker_callback4, button='Middle')
