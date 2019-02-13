import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches
from snowav import database
from snowav.database.tables import Basins
from snowav.plotting.plotlims import plotlims as plotlims
import pandas as pd
from matplotlib.ticker import FormatStrFormatter


def subbasins(snow):


    lims = plotlims(snow.basin, snow.plotorder)

    demcopy = copy.copy(snow.dem)
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0
    demcopy[ixo] = np.nan

    # filter out high and low hyp
    map = plt.cm.terrain
    map.set_bad('white')

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    plt.close(31)
    f = plt.figure(num=31, figsize=(6,6), dpi=snow.dpi)
    ax = plt.gca()

    h = ax.imshow(demcopy, cmap=map, alpha=0.45)

    # Basin boundaries
    for name in snow.masks:
        ax.contour(snow.masks[name]['mask'], cmap='Greys',linewidths=1)

    # If these attributes are included in config, put labels on ax
    if (snow.annot_x is not None) and (snow.annot_y is not None):

        for i, name in enumerate(lims.sumorder):
            if name in ['Willow Creek','South Fork',
                        'Upper South Fork','Lower South Fork',
                        'Middle South Fork', 'Kaweah River',
                        'Lake Kaweah', 'Dinkey Creek']:
                if len(name.split(' ')) == 2:
                    n = '{}\n{}'.format(name.split(' ')[0],name.split(' ')[1])
                else:
                    n = '{}\n{} {}'.format(name.split(' ')[0],name.split(' ')[1],name.split(' ')[2])
            else:
                n = name

            ax.annotate(n, xy=(snow.annot_x[i], snow.annot_y[i]), size = 12)

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])

    plt.tight_layout()

    plt.savefig('{}subbasins_{}.png'.format(snow.figs_path,snow.name_append))
