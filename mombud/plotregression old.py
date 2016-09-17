# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 00:16:45 2016

@author: sweel_Rafelski
"""


def getrval(df, x, y, labeldic):
    """
    return a subset DataFrame and R^2 values for columns x, y in original df
    """
    df = df[[x, y, 'date']].reset_index(drop=True)
    df.rename(columns=labeldic, inplace=True)
    pear = (df.groupby('date')
            .corr().xs(labeldic[y], level=1)
            ['raw'].to_dict())
    r_sqr = {key: value**2 for key, value in pear.iteritems()}

    return df, r_sqr


def plotRegr(**kwargs):
    """
    plot regression coeff of Δψ raw vs scaled
    """
    df = kwargs.get('data')
    datadir = kwargs.get('savefolder')
    save = kwargs.get('save', False)

    labeldic = {'DY_abs_cell_mean': 'raw',
                'DY_cell_mean': 'scaled cell mean',
                'DY_median_mom': 'scaled mom mean',
                'DY_median_bud': 'scaled bud mean'}

    with sns.plotting_context('talk', font_scale=1.):
        for p in ['DY_cell_mean', 'DY_median_mom', 'DY_median_bud']:
            a, r2 = getrval(df, 'DY_abs_cell_mean', p, labeldic)
            x, y, _ = a.columns
            (sns.lmplot(x, y,
                        fit_reg=False,
                        legend_out=False,
                        data=a,
                        hue='date',
                        size=8,
                        aspect=1.5)
             .set(xlim=(0, 2000), xlabel='raw GFP'))
            if save:
                plt.savefig(op.join(datadir,
                                    "{} vs raw.png".format(labeldic[p])))

            g = (sns.lmplot(x, y,
                            fit_reg=False,
                            data=a,
                            col='date',
                            hue='date',
                            col_wrap=3)
                 .set(xlim=(0, 2000), ylim=(0., 1.), xlabel='raw GFP'))

            labelRsqr(g, r2)

            if save:
                plt.savefig(op.join(datadir,
                                    "{} vs raw by date.png".format(
                                        labeldic[p])))


def ranger(lst, mult):
    """
    return array for y-ticks or x-ticks with multiples of `mult`
    from limits in `lst`
    """
    yt = [i - i % mult for i in lst]
    return np.arange(yt[0], yt[1], mult)


def plotDyAxisDist(dfmom, dfbud, **kwargs):
    """
    Progression of Δψ as move along the bud axis
    """

    bigbinsbud = pd.melt(dfbud,
                         id_vars=['media', 'binvol'],
                         var_name='bud axis position',
                         value_name=r'$\Delta\Psi$ scaled gradient',
                         value_vars=binsaxis.tolist())
    bigbinsbud = bigbinsbud.dropna()

    # with facetting by budvol
    with sns.plotting_context('talk', font_scale=.9):
        ylims2 = (bigbinsbud[r'$\Delta\Psi$ scaled gradient']
                  .quantile([0.025, 0.975]).tolist())

        if (ylims[1] - ylims[0]) < 2:
            mult = 0.25
        else:
            mult = np.ceil((ylims2[1] - ylims2[0]) / 6)
        yt = ranger(ylims2, mult)

        m0 = sns.FacetGrid(bigbinsbud,
                           row='media',
                           col="binvol",
                           hue='media',
                           row_order=COL_ODR,
                           col_order=binsvolbud[1:])

        m0 = (m0.map(sns.pointplot,
                     'bud axis position',
                     r'$\Delta\Psi$ scaled gradient')
              .set(yticks=yt, ylim=tuple(ylims2)))
        labelBudVol(m0, Nbud)

        if save:
            m0.savefig(op.join(datadir, 'bud_cell_dy_facetted.png'))

        # =========================================================================
        # Plotting routines
        # =========================================================================
        print "List abbr. for plotting function names\n"
        abbrv_names, funcdict = getFuncList(vp, 'plot')
        for f in abbrv_names:
            print "{key}: {func}".format(key=f, func=funcdict[f].__name__),
            print "{docs}".format(docs=funcdict[f].__doc__)

        while True:
            try:
                invar = raw_input('Please enter abbrev. two letter name '
                                  'of functions to plot,\n'
                                  'or "a" for all plots, '
                                  'or "q" to quit (without quotes): ')

                if invar == 'q':
                    print 'Quitting!'
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

        return 0
