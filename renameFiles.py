import os
import os.path as op
import glob
import string
import shutil as sh
import fnmatch as fn
import errno


def copy(src, dest):
    try:
        sh.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            sh.copy(src, dest)
        else:
            print('Directory not copied. Error: %s' % e)

parDir = os.path.dirname(os.getcwd())
b = glob.glob(parDir+'\*[0-9][0-9]*')  # directory names
A = glob.glob('N*.vtk')
copyTo = b[0]
toMove = ['YPD_', 'YPE_', 'YPR_', 'YPL_', 'YPL-']
for g, h in enumerate(toMove):
    for temp in A:
        if (temp[5:9] == h):
            sh.move(temp, b[g])

# =============================================================================
#                  filter for a filename and get the path
# =============================================================================
temp = []
for root, dirs, files in os.walk(os.getcwd()):
    for f in files:
        if fn.fnmatch(f, 'zMax*'):
            temp.append(
                os.path.join(root, f))

# get the root folders first  before getting the list of stack files
temp = []
for root, dirs, files in os.walk(os.getcwd()):
    for f in dirs:
        temp.append(
            os.path.join(root, f))
temp = temp[:3]

#  after that filter and move to new directory
for i in temp:
    for root, dirs, files in os.walk(i):
        try:
            for i in files:
                if fn.fnmatch(i, '*RFPs*tif'):
                    sh.move(
                        os.path.join(
                            root, i),
                        os.path.join(
                            root, 'rfpStacks', i))
        except:
            IOError

# =============================================================================
#       find skel files and copy resampled file with
#       corresponding name to one directory
# =============================================================================
for root, dirs, files in os.walk(os.getcwd()):
    for f in files:
        if fn.fnmatch(f, '*skel*vtk'):

            temp2 = f.replace("RFP", "GFP")

            resampGFP = string.join(
                (temp2[:-13], "resampled.vtk"), sep='_')

            sh.copy(resampGFP,
                        os.path.join(root, "toKeep"))  # dir to store all files

# =============================================================================
# MOVEs IMAGES ONE FOLDER UP!!! note the inplace mod of dirs[:] and the bottom
#up traversal in the first walk level!!
#A similar result could have been done with storing the directories in a temp
#list first and then walking through each directory in the list
# =============================================================================

# ****** RUN THESE FILES IN ORDER FIRST ***************************
# 1) move up one folder from 'Pos'
for root, dirs, _ in os.walk(os.getcwd(), topdown=False):
    dirs[:]=[d for d in dirs if fn.fnmatch(d, "*Pos*")]  # inplace list modfhghfghgh
    if dirs:
        for r, d, files in os.walk(op.join(root, dirs[0])):
            for f in files:
                sh.move(op.join(root, dirs[0], f),
                        op.join(root, f))

# 2) repad the FOLDER names middle index to 00[0-9] format
for dirs in os.listdir(os.getcwd()):
    olddir = os.path.abspath(dirs)
    oldstr = os.path.basename(olddir)
    oldstrL = oldstr.split("_")
    newstr = "_".join((oldstrL[0], oldstrL[1].zfill(3)))
    newdir =  olddir.replace(oldstr, newstr)
    os.rename(olddir, newdir)

# 3) rename the files to switch labels
for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fn.fnmatchcase(i, '*zQUAD GFP*tif'):
            i2 = op.join(root,i)
            newpath = i2.replace(" GFP", "_rfp")
            os.rename(op.join(root, i),
                      op.join(newpath))

# optional delete if mess up from #3
for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fn.fnmatchcase(i, '*_rfp_*tif'):
            print op.join(root,i)
            os.remove(op.join(root, i))
