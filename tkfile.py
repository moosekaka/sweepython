# -*- coding: utf-8 -*-
"""
Created on Thu May 26 13:46:43 2016

@author: sweel_Rafelski
"""

import Tkinter as tk
from tkFileDialog import *

root = tk.Tk()
root.withdraw()
file_path = askopenfilename(title ='Select a file!')