# -*- coding: utf-8 -*-
"""
Created on Fri Feb 03 10:52:39 2017

"""

import Tkinter, tkFileDialog

class SelectDirClient(Tkinter.Frame):

  def __init__(self, root, **kwargs):
    Tkinter.Frame.__init__(self, root)
    self.root = root
    self.path = None
    # defining options for opening a directory
    self.dir_opt = options = {}
    options['initialdir'] = kwargs.pop('initialdir', '.')
    options['parent'] = root
    options['title'] = kwargs.pop('title', 'Select Folder')

  def askdirectory(self):

    """Returns a selected directoryname."""
    self.path = tkFileDialog.askdirectory(**self.dir_opt)
    self.root.destroy()
#    return path

if __name__=='__main__':
  root = Tkinter.Tk()
  root.withdraw()
  gui = SelectDirClient(root)
  gui.pack()
  root.after(500, gui.askdirectory)
  root.mainloop()

