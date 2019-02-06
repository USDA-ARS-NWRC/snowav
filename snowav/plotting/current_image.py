
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import copy
import cmocean
import matplotlib.patches as mpatches

def current_image(snow):
    '''
    This plots images of SWE and cold content at snow.end_date.

    '''

    state = copy.deepcopy(snow.outputs['swe_z'][snow.ixe])
    state = np.multiply(state,snow.depth_factor)
    cold = copy.deepcopy(snow.outputs['coldcont'][snow.ixe])
    cold = np.multiply(cold,0.000001)

    qMin,qMax = np.nanpercentile(state,[0,99.8])
    clims = (qMin,qMax)
    clims2 = (-5,0)
    pmask = snow.masks[snow.plotorder[0]]['mask']
    ixo = pmask == 0

    # Prepare no-snow and outside of the basin for the colormaps
    ixz = state == 0
    state[ixz] = -1
    cold[ixz] = 1

    # Colormap for snow.state
    colorsbad = plt.cm.Set2_r(np.linspace(0., 1, 1))
    colors1 = cmocean.cm.haline_r(np.linspace(0., 1, 254))
    colors = np.vstack((colorsbad,colors1))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    state[ixo] = np.nan
    mymap.set_bad('white',1.)

    mymap1 = plt.cm.Spectral_r
    cold[ixo] = np.nan
    mymap1.set_bad('white')
    mymap1.set_over('lightgrey',1)

    sns.set_style('dark')
    sns.set_context("notebook")

    plt.close(1)
    fig,(ax,ax1) = plt.subplots(num=1, figsize=snow.figsize,
                                facecolor = 'white', dpi=snow.dpi,
                                nrows = 1, ncols = 2)
    h = ax.imshow(state, cmap = mymap, clim=clims)
    h1 = ax1.imshow(cold, clim=clims2, cmap = mymap1)

    # If these attributes are included in config, put labels on ax
    if (snow.annot_x is not None) and (snow.annot_y is not None):
        for i, name in enumerate(snow.plotorder[1::]):
            if name in ['Willow Creek','South Fork',
                        'Upper South Fork','Lower South Fork',
                        'Middle South Fork']:
                if len(name.split(' ')) == 2:
                    n = '{}\n{}'.format(name.split(' ')[0],name.split(' ')[1])
                else:
                    n = '{}\n{} {}'.format(name.split(' ')[0],name.split(' ')[1],name.split(' ')[2])
            else:
                n = name

            ax.annotate(n, xy=(snow.annot_x[i], snow.annot_y[i]), size = 8)

    if snow.basin == 'LAKES':
        ax.set_xlim(snow.imgx)
        ax.set_ylim(snow.imgy)
        ax1.set_xlim(snow.imgx)
        ax1.set_ylim(snow.imgy)

    # Basin boundaries
    for name in snow.masks:
        ax.contour(snow.masks[name]['mask'],cmap = "Greys",linewidths = 1)
        ax1.contour(snow.masks[name]['mask'],cmap = "Greys",linewidths = 1)

    # if snow.basin == 'SJ':
    #     fix1 = np.arange(1275,1377)
    #     fix2 = np.arange(1555,1618)
    #     ax.plot(fix1*0,fix1,'k')
    #     ax.plot(fix2*0,fix2,'k')
    #     ax1.plot(fix1*0,fix1,'k')
    #     ax1.plot(fix2*0,fix2,'k')

    # Do pretty stuff for the left plot
    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cbar.set_label('[%s]'%(snow.depthlbl))

    # Do pretty stuff for the right plot

    h1.axes.get_xaxis().set_ticks([])
    h1.axes.get_yaxis().set_ticks([])
    h1.axes.set_title('Cold Content \n %s'%(snow.report_date.date().strftime("%Y-%-m-%-d")))
    divider = make_axes_locatable(ax1)
    cax2 = divider.append_axes("right", size="5%", pad=0.2)
    cbar1 = plt.colorbar(h1, cax = cax2)
    cbar1.set_label('Cold Content [MJ/$m^3$]')
    cbar1.ax.tick_params()

    h.axes.set_title('SWE \n %s'%(snow.report_date.date().strftime("%Y-%-m-%-d")))
    fig.subplots_adjust(top=0.95,bottom=0.05,
                        right = 0.92, left = 0.05, wspace = 0.12)

    if snow.basin == 'LAKES':
        plt.tight_layout()

    patches = [mpatches.Patch(color='grey', label='snow free')]

    # Make legend/box defaults and adjust as needed
    pbbx = 0.3
    pbby = 0.05

    if snow.basin in ['KAWEAH', 'RCEW']:
        pbbx = 0.1

    # snow-free
    ax.legend(handles=patches, bbox_to_anchor=(pbbx, pbby),
              loc=2, borderaxespad=0. )

    snow._logger.info('saving figure to %sresults_%s.png'%(snow.figs_path,snow.name_append))
    plt.savefig('%sresults_%s.png'%(snow.figs_path,snow.name_append))
