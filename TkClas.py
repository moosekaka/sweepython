# -*- coding: utf-8 -*-
"""
Created on Fri Feb 03 10:52:39 2017

"""

import Tkinter, Tkconstants, tkFileDialog

class SelectDirClient(Tkinter.Frame):

  def __init__(self, root, **kwargs):

    Tkinter.Frame.__init__(self, root)

    # options for buttons
    button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}

    # define buttons
    Tkinter.Button(self, text='askdirectory',
                   command=self.askdirectory).pack(**button_opt)

    # defining options for opening a directory
    self.dir_opt = options = {}
    options['initialdir'] = kwargs.pop('initialdir', '.')
    options['parent'] = root
    options['title'] = kwargs.pop('title', 'Select Folder')

  def askdirectory(self):

    """Returns a selected directoryname."""

    return tkFileDialog.askdirectory(**self.dir_opt)

if __name__=='__main__':
  root = Tkinter.Tk()
  SelectDirClient(root).pack()
  root.mainloop()

