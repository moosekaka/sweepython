# -*- coding: utf-8 -*-
"""
Main module to analyze mom bud asymmetry
"""
import os
import os.path as op
import cPickle as pickle
import inspect
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mombud.functions import vtk_mbfuncs as vf
from mombud.functions import vtk_mbplots as vp
from wrappers import UsageError, swalk, ddwalk

# pylint: disable=C0103


def getFuncList(module, splitword):
    """
    generate a dispatch function hash with two letter abbreviation of the
    functions name as key
    """
    fun_list = [f for f in inspect.getmembers(module, inspect.isfunction)
                if f[0].startswith(splitword)]

    fun_dict = {}
    for tup in fun_list:
        abbr = tup[0].split('plot')[1][:2].lower()
        fun_dict[abbr] = tup[1]
    return sorted(fun_dict.keys()), fun_dict


def getData():
    """
    Get input data from specified work dirs

    Returns
    -------
    filekeys_f, datadir : dict(Str), Str
        filepaths of inputs and folder for data files
    dfmb : DataFrame
        cell volume data
    """

    datadir = op.join(os.getcwd(), 'mutants', 'transformedData', 'filtered')

    # old data
    datadir_old = op.join(os.getcwd(), 'data', 'transformedData')

    # DataFrames for new and old cell picked point
    try:
        with open(op.join(datadir, 'mombudtrans_new.pkl'), 'rb') as inpt:
            dfmb = pickle.load(inpt)  # cols base, neck, tip, media, bud, mom

        with open(op.join(datadir_old, 'mombudtrans.pkl'), 'rb') as inpt:
            dfmb_o = pickle.load(inpt)

        # reject candidates
        rejectfold = op.join(datadir, os.pardir, 'reject')
        reject = swalk(rejectfold, '*png', stop=-4)

    # VTK files for new and old data
        filext = "*vtk"
        vtkF = ddwalk(datadir, filext, stop=-4)
        vtkF_old = swalk(datadir_old, filext, stop=-4)

    except UsageError as e:
        print e.args[0]
        raise

    # file paths for VTKs
    filekeys_old = {item: vtkF_old[item] for item
                    in sorted(vtkF_old.keys()) if
                    item.split('_')[0] != 'YPD' and
                    item not in reject}

    filekeys = {item: vtkF[media][item] for media
                in sorted(vtkF.keys()) for item
                in sorted(vtkF[media].keys())}

    filekeys_f = {f: filekeys[f] for f in filekeys if f not in reject}
    filekeys_f.update(filekeys_old)  # add YPE to dict
    # dataframe of neck, mom and bud tip positions, bud and mom volumes
    dfmb = dfmb.append(dfmb_o[dfmb_o.media != 'YPD'])
    return filekeys_f, dfmb, datadir


def _concatDF(vtkdf):
    keys = sorted(vtkdf.keys())
    dic = defaultdict(dict)
    for k in keys:
        # cell is ref/view (not deep copy) of vtkdf[k]['df'], changes to cell
        # results in changes to vtkdf[k]['df'], DataFrames are mutable
        cell = vtkdf[k]['df']
        # update with whole cell stat. data
        dic['cell'][k] = vtkdf[k]['celldata']
        # set index to cellname so that concatenate is possible
        cell['name'] = k
        cell.set_index('name', inplace=True)
    # Concat of ALL cell dfs into one giant DataFrame
    df_concat = pd.concat([vtkdf[k]['df'] for k in keys])
    df_concat.reset_index(inplace=True)
    return df_concat, dic


def _scaleDY(df):
    """
    scaling for group date variations in Δψ
    """
    grd = df.groupby('date')
    # normalize by date mean DYunscl
    df['DYun_f'] = (grd['DY_unscl']
                    .transform(lambda x: (x - x.mean()) / x.std()))
    df['DYun_f2'] = (grd['DY_unscl']
                     .transform(lambda x: x - x.mean()))
    df['DYun_f3'] = (grd['DY_unscl']
                     .transform(lambda x: (x - x.min()) / (x.max() - x.min())))


def _aggDY(df):
    gr = df.groupby('name')
    labels = gr.first()[['date', 'media']]
    # groupby mom/buds , get agg. stats for Δψ
    df_agg_mb = (df.groupby(['name', 'type'])
                 [['DY', 'DY_abs', 'DY_unscl',
                   'DYun_f', 'DYun_f2', 'DYun_f3']]
                 .agg([np.mean, np.median]).unstack())
    df_agg_mb.columns = ['_'.join(c) for c in df_agg_mb.columns.values]

    df_agg_cell = (gr[['DY', 'DY_abs', 'DY_unscl',
                       'DYun_f', 'DYun_f2', 'DYun_f3']].agg('mean'))
    df_agg_cell.rename(columns=lambda x: x + '_cell_mean', inplace=True)

    return pd.concat([df_agg_mb, df_agg_cell, labels], axis=1)


def _mombudDF(df, dic, dy_type='DY'):
    """
    groupby bins of ind cell position
    """
    gr = df.groupby(['name', 'type', 'ind_cell_binpos'])
    dfbinned = (gr[dy_type].mean()
                .unstack(level='ind_cell_binpos'))
    dfbinned.columns = dfbinned.columns.astype('float')

    # scale by whole cell mean Δψ
    df = dic['dfcell']['DY_cell_mean']
    for i in ['dfbud', 'dfmom']:
        dic[i] = dfbinned.xs(i[2:], level='type')
        dic[i] = dic[i].div(df, axis=0)  # scaling by whole cell mean Δψ


def _update_dfMB(key, df_all, df_ind, **kwargs):
    """
    concatenate type and vol. data from df_all df to ind. mom/bud df_ind
    """
    bins = kwargs['binsvol%s' % key]
    df_ind = df_ind.assign(media=df_all.loc[:, 'media'],
                           date=df_all.loc[:, 'date'])
    df_ind['%svol' % key] = df_all['%svol' % key]
    df_ind['binvol'] = pd.cut(df_ind['%svol' % key],
                              bins=bins, labels=bins[1:])
    return df_ind


def _subsetDF(df, dic, keylist=None):
    """
    return a subset, eg. YPE data from full df dataset
    """
    if keylist is None:
        keylist = ['YPE', 'WT']
    subset = df[df['media'].isin(keylist)]
    subset = subset.reset_index(drop=True)
    cntlab = keylist[0].lower()
    dic['counts_%s' % cntlab] = subset.groupby('date').size().to_dict()
    dic['data_%s' % cntlab] = subset


def _filterMask(df):
    """
    Filter conditions, reject large cells
    """
    filt_large_ypd = (df.momvol > 100) & (df.media != 'YPD')
    filt_type = (((df.media == 'YPE') & (df.date == '052315')) |
                 ((df.media == 'WT') & (df.date == '032716')))
    maskdic = {'large_ypd': ~(filt_large_ypd)}
    maskdic['large_ypd_hilo_ype'] = ~(filt_large_ypd) & ~(filt_type)
    return maskdic


def process_ind_df(vtkdf, mbax=None, cellax=None, **kwargs):
    """
    bin Δψ distrbution along cellaxis for each ind. cell DataFrame and append
    to the respective DataFrames

    Parameters
    ----------
    vtkdf : DataFrame
        inndividual DataFrame inputs to be concatenated
    mbax : np.array
        bin cell position for ind. mom/bud cell
    cellax : np. array
        bin cell position for whole cell
    Returns
    -------
    dicout : dict
        dictionary of DataFrames for mom bud analyses
    """

    if mbax is None:
        raise UsageError('please specify kwargs for "mbax" '
                         '(mom bud axis bins)')

    if cellax is None:
        raise UsageError('please specify kwargs for "cellax" '
                         '(whole cell axis bins)')

    # concat individual cell DFs and get individual cell data dic_in
    dfc, dic_in = _concatDF(vtkdf)

    # bin the dataframe according to individual (mom/bud) axis
    groups = dfc.groupby('name')
    dfc['ind_cell_binpos'] = (groups['ind_cell_axis']
                              .apply(pd.cut, bins=mbax,
                                     labels=mbax[1:]))

    # bin the dataframe according to individual entire cell axis
    dfc['whole_cell_binpos'] = (groups['whole_cell_axis']
                                .apply(pd.cut, bins=cellax,
                                       labels=cellax[1:]))

    # get date and mediatype str labels
    split = dfc['name'].str.split('_')
    dfc['media'] = [x[0] for x in split]
    dfc['date'] = [x[1].replace('c', '0') if x[1].startswith('c') else x[1]
                   for x in split]

    # Calc. scaling factor for raw GFP daily variations
    _scaleDY(dfc)
    # aggregated mean Δψ by date and by cell
    dfc_agg = _aggDY(dfc)

    # DataFrame for agg. mean data of ind. cells
    dicout = defaultdict(dict)
    dicout['dfcell'] = pd.DataFrame.from_dict(dic_in['cell'], orient='index')
    dicout['dfcell'] = pd.concat([dicout['dfcell'], dfc_agg], axis=1)
    dicout['concat'] = dfc
    _mombudDF(dfc, dicout, dy_type='DY')
    return dicout


def postprocess_df(**kwargs):
    """
    Set population level data ,update parameters dict for plotting and filter
    unwanted data

    Returns
    -------
    outputdic : dict
        dictionary of DataFrames for mom bud analyses
    """

    kwargs['filekeys_f'], kwargs['dfmb'], kwargs['savefolder'] = getData()
    ind_cell_df = vf.gen_data(dfvoldata=kwargs['dfmb'],
                              fkeys=kwargs['filekeys_f'], **kwargs)
    Dout = process_ind_df(ind_cell_df, **kwargs)
    cellall = Dout.pop('dfcell')
    cellposmom = Dout.pop('dfmom')
    cellposbud = Dout.pop('dfbud')

    # add cell volume data from cell tracing data
    cellall['budvol'] = kwargs['dfmb'].bud
    cellall['momvol'] = kwargs['dfmb'].mom

    # frac -> ratio of mom/bud Δψ
    cellall = (cellall
               .assign(frac=cellall.loc[:, 'DY_median_bud'] /
                       cellall.loc[:, 'DY_median_mom']))

    #  normalize budvol_q90 -> largest cells (90th percentile)
    cellall = (cellall
               .assign(budvol_q90=cellall['media']
                       .map(cellall.groupby('media')['budvol']
                            .quantile(.90))))

    cellall['budvolratio'] = (cellall['budvol']
                              .div(cellall['budvol_q90'], axis=0))

    # Output dict. for cellall, Δψ binned by mom and bud ind. cells
    outputdic = {'data': cellall}  # for all cells
    outputdic['dfmom'] = _update_dfMB('mom', cellall, cellposmom, **kwargs)
    outputdic['dfbud'] = _update_dfMB('bud', cellall, cellposbud, **kwargs)

    # Bins used for plotting budding progression
    binsaxisbig = kwargs['cellax']  # 2. cat. for cells > 90th percentile
    binsaxisbig = np.r_[binsaxisbig, [2.]]
    cellall['bin_budprog'] = pd.cut(cellall['budvolratio'],
                                    bins=binsaxisbig, labels=binsaxisbig[1:])
    cellall['binbudvol'] = outputdic['dfbud']['binvol']

    # filter out criteria
    filtout = _filterMask(cellall)
    for i in ['data', 'dfmom', 'dfbud']:
        outputdic[i] = outputdic[i][filtout['large_ypd']]

    # get counts for each type
    outputdic['counts'] = (outputdic['data']
                           .groupby('media')
                           .size().to_dict())
    outputdic['counts_buds'] = (outputdic['data']
                                .groupby(['media', 'binbudvol'])
                                .size())
    outputdic['counts_date'] = (outputdic['data']
                                .groupby('date')
                                .size().to_dict())

    # subset for YPE
    _subsetDF(cellall, outputdic)

    # output as dict.
    outputdic.update(kwargs)
    outputdic.update(Dout)  # any leftover vars in Dout are returned
    return outputdic


def main(**kwargs):
    """
    Main

    kwargs
    ------
    plotlist : List
        plotting function names
    """

    os.chdir(op.expanduser(os.sep.join(
        ('~', 'Documents', 'Github', 'sweepython', 'WorkingData'))))

    def_args = {'regen': False,
                'save': False,  # toggle to save plots
                'inpdatpath': 'celldata.pkl',
                'mbax': np.linspace(0., 1., 6),  # pos. along mom/bud cell
#                'cellax': np.linspace(0, 1., 11),  # position along whole cell
                'binsvolbud': np.linspace(0, 40, 5),  # vol binning for bud
                'binsvolmom': np.array([0, 30, 40, 80.]),
                'gfp_plot_vars': ['DY_abs_mean_mom',
                                  'DY_abs_mean_bud',
                                  'DY_abs_cell_mean'],  # plotGFP variables
                'viol_plot_vars': ['DY_median_mom',
                                   'DY_median_bud',
                                   'DY_abs_cell_mean'],  # plotviol variables
                'COL_ODR': ['MFB1', 'NUM1', 'YPT11',
                            'WT', 'YPE', 'YPL', 'YPR']}

    def_args.update(kwargs)  # override default args with user kwargs, if any
    outputargs = postprocess_df(**def_args)  # call getdata(), process_ind_df()

    # =========================================================================
    # Plotting routines
    # =========================================================================
    print "List abbr. for plotting function names\n"
    abbrv_names, funcdict = getFuncList(vp, 'plot')
    for f in abbrv_names:
        print "{key}: {func}".format(key=f, func=funcdict[f].__name__),
        print "{docs}".format(docs=funcdict[f].__doc__)

    plot_switch = kwargs.get('plot_switch', True)
    while plot_switch:
        try:
            invar = raw_input('Please enter abbrev. two letter name '
                              'of functions to plot,\n'
                              'or "a" for all plots, '
                              'or "q" to quit (without quotes): ')

            if invar == 'q':
                break

            if invar == 'a':
                _ = [funcdict[f](**outputargs) for f in abbrv_names]
                print "Finished plotting {} functions!".format(len(_))
                break

            funcdict[invar](**outputargs)
            print "{} has been executed!".format(funcdict[invar].__name__)

        except KeyError as e:
            print "{} is not in functions list".format(e)
            continue

# _____________________________________________________________________________
if __name__ == '__main__':
    plt.close('all')
    main(regen=False, plot_switch=True, save=True)
