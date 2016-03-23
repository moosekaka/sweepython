# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 18:07:36 2015

@author: sweel
"""

figure.scene.disable_render = False  # also something new i've learnt!


def picker_callback(picker):
    """ Picker callback: this get called when on pick events.
    """
    global mypts_selection, myplane
    glyph_points = mypts.glyph.glyph_source.glyph_source.output.points.to_array()
    point_id = picker.point_id / glyph_points.shape[0]
    # If the no points have been selected, we have '-1'
    if point_id != -1:
        # Retrieve the coordinnates coorresponding to that data
        # point
        x, y, z = xtr[point_id, 0], xtr[point_id, 1], xtr[point_id, 2]
        # Move the outline to the data point.
        outline.bounds = (x - 0.1, x + 0.1,
                          y - 0.1, y + 0.1,
                          z - 0.1, z + 0.1)
        pt = xtr[point_id, 0:3]
        pts = tree.query_ball_point(pt, r=1.)

        figure.scene.disable_render = True
        # draw the nearest few points to the selected one as red dots
        if mypts_selection is not None:
            mypts_selection.parent.parent.remove()
        mypts_selection = mlab.points3d(xtr[pts, 0], xtr[pts, 1], xtr[pts, 2], scale_factor=.15, color=(
            1, 0, 0), mode='sphere', scale_mode='none', colormap='prism')  # ,mask_points=10)

        # fit a plane to the point neighbourhoud
        B, normd = fitplane(xtr[pts, 0:3].T)
        # and draw that plane
        if myplane is not None:
            myplane.remove()
        myplane = mlab.pipeline.builtin_surface()
        myplane.source = 'plane'
        myplane.data_source.normal = B[0:3]
        myplane.data_source.center = Plane(P=Point(B[0], B[1], B[2]), D=-B[3]).projection(
            Point(xtr[point_id, 0], xtr[point_id, 1], xtr[point_id, 2])).asarray()
        planesurface = mlab.pipeline.surface(myplane)
        figure.scene.disable_render = False

# set the picker callback function and mouse tolerance
picker = figure.on_mouse_pick(picker_callback)
picker.tolerance = 0.01

# start the magic
mlab.show()
