import os
import os.path as op
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


def mkdir_exist(path):
    try:
        os.makedirs(path)
    except OSError as error:
        if error.errno != errno.EEXIST:
            raise


def renamecopy(f, r, word):
    """
    macro to change filename label to date and return path of file
    """
    if fn.fnmatch(f, word):
        filename = f.partition("_")
        pth = op.join(op.dirname(r), filename[0])
        date = r.rsplit(os.sep)[-1]
        newname = "_".join((date, filename[-1]))
        return (pth, newname)


# =============================================================================
# MOVEs IMAGES ONE FOLDER UP!!! note the inplace mod of dirs[:] and the bottom
# up traversal in the first walk level!!
# A similar result could have been done with storing the directories in a temp
# list first and then walking through each directory in the list
# =============================================================================

# ****** RUN THESE FILES IN ORDER FIRST ***************************
# 1) move up one folder from 'Pos'
for root, dirs, _ in os.walk(os.getcwd(), topdown=False):
    dirs[:] = [d for d in dirs if fn.fnmatch(
        d, "*Pos*")]  # inplace list modfhghfghgh
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
    newdir = olddir.replace(oldstr, newstr)
    os.rename(olddir, newdir)

# 3) rename the files to switch labels
for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fn.fnmatchcase(i, '*zQUAD GFP*tif'):
            i2 = op.join(root, i)
            newpath = i2.replace(" GFP", "_rfp")
            os.rename(op.join(root, i),
                      op.join(newpath))

# optional delete if mess up from #3
for root, dirs, files in os.walk(os.getcwd()):
    for i in files:
        if fn.fnmatchcase(i, '*_rfp_*tif'):
            print op.join(root, i)
            os.remove(op.join(root, i))


# =============================================================================
# MUTANT CELLS Dataset (YPE)
# =============================================================================
# copy vtk files into app. folders

for root, dirs, files in os.walk(os.getcwd(), topdown=False):
    # search for foldernames starting with two numbers
    dirs[:] = [op.join(root, d) for d in dirs if fn.fnmatch(d, "*[0-9][0-9]*")]

for dircell in dirs:
    for root, dirs, files in os.walk(dircell):
        for fil in files:
            P = renamecopy(fil, root, '*RFP*skeleton*')  # change
            if P:
                mkdir_exist(P[0])
                sh.copy(op.join(root, fil), op.join(P[0], P[1]))


# for root, dirs, files in os.walk(os.getcwd()):
#    for f in files:
#        if fn.fnmatchcase(f, '*surface*'):
#            old = op.join(root, f)
#            new = old + ".vtk"
#            os.rename(old, new)

# =============================================================================
# rename BF stacks to be same name as RFP for mom bud tracing macros
# =============================================================================
for root, dirs, files in os.walk(os.getcwd(), topdown=False):
    dirs[:] = [op.join(root, d) for d in dirs if fn.fnmatch(d, "[BGR]F*")]

BF = [files for _, _, files in os.walk(dirs[0])][0]
RFP = [files for _, _, files in os.walk(dirs[2])][0]

newname = {}
for key, val in zip(BF, RFP):
    newname[key] = val

for bf in BF:
    oldpath = op.join(root, 'BF', bf)
    newpath = op.join(root, 'BF', newname[bf])
    os.rename(oldpath, newpath)
