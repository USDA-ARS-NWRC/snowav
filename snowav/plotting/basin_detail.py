import math
import copy
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from spatialnc import ipw
from shutil import copyfile
import os
import pandas as pd


def basin_detail(self):

    demcopy = copy.copy(self.dem)
    pmask = self.masks[self.plotorder[0]]['mask']
    ixo = pmask == 0
    demcopy[ixo] = np.nan

    emin = int(math.ceil(np.nanmin(demcopy)/10.0))*10
    emax = int(math.ceil(np.nanmax(demcopy)/10.0))*10
    step = 50 # 50
    edges = np.arange(emin,emax+step,step)
    ixd = np.digitize(demcopy,edges)
    hypsom = pd.DataFrame(index = edges, columns = self.masks.keys())

    for name in self.masks:
        elevbin = np.multiply(ixd,self.masks[name]['mask'])

        for n in np.arange(0,len(edges)):
            ind        = elevbin == n
            hypsom.loc[edges[n],name] = (np.nansum(np.nansum(ind))/np.nansum(np.nansum(self.masks[name]['mask'])))*100

    for name in self.masks:
        smin = np.nanmin(self.dem[self.dem*self.masks[name]['mask'] > 0])
        smax = np.nanmax(self.dem[self.dem*self.masks[name]['mask'] > 0])
        self._logger.info(name,str(int(smin),str(int(smax))))

    # filter out high and low hyp
    map = plt.cm.terrain
    map.set_bad('white')

    sns.set_style('darkgrid')
    sns.set_context("notebook")

    clrs = copy.copy(self.barcolors)
    clrs.insert(0,'k')

    plt.close(10)
    fig,(ax,ax1) = plt.subplots(num=10, figsize=self.figsize, dpi=self.dpi, nrows = 1, ncols = 2)
    h = ax.imshow(demcopy, cmap=map)

    # Basin boundaries
    for name in self.masks:
        ax.contour(self.masks[name]['mask'],cmap = "Greys",linewidths = 1)

    for iters,name in enumerate(self.plotorder):
        hyp = copy.deepcopy(hypsom[name].cumsum()-hypsom[name].iloc[0])
        imin = hyp < 0.1
        hyp[imin] = np.nan
        imax = hyp > 99.9
        hyp[imax] = np.nan
        ax1.plot(hyp,range(0,len(edges)),color = clrs[iters],label = name)

    ax1.invert_xaxis()
    ax1.legend()
    ax1.set_xlim((102,-2))
    ax1.set_ylim((-1,len(hypsom)))

    ax1.set_ylabel('elevation [%s]'%(self.elevlbl))
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")

    labels = ax1.get_xticklabels()
    ax1.set_xticklabels(labels[::-1])

    xts = ax1.get_yticks()
    edges_lbl = []
    for i in xts[0:len(xts)-1]:
        edges_lbl.append(str(int(edges[int(i)])))

    ax1.set_yticklabels(str(i) for i in edges_lbl)

    ax1.set_xlabel(r'Basin area above elevation [%]')
    ax.set_title(self.plotorder[0])

    h.axes.get_xaxis().set_ticks([])
    h.axes.get_yaxis().set_ticks([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("left", size="5%", pad=0.2)
    cbar = plt.colorbar(h, cax = cax)
    cax.yaxis.set_ticks_position('left')
    cax.yaxis.set_label_position('left')
    cax.set_ylabel('elevation [m]')

    plt.subplots_adjust(top=0.88)
    plt.tight_layout()

    if self.basin == 'SJ':
        ax.text(0.3,0.88,'Main',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
        ax.text(0.48,0.055,'Jose Creek',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
        ax.text(0.95,0.65,'South\nFork',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
        ax.text(0.05,0.55,'Willow\nCreek',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
    if self.basin == 'BRB':
        ax.text(0.2,0.74,'Mores Creek',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
        ax.text(0.48,0.77,'Twin\nSprings',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
        ax.text(0.88,0.6,'Featherville',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
    if self.basin == 'Tuol':
        ax.text(0.75,0.68,'Tuolumne',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
        ax.text(0.2,0.84,'Cherry',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)
        ax.text(0.14,0.35,'Eleanor',horizontalalignment='center',transform=ax.transAxes,fontsize = 10)

    self._logger.info('saving figure to %shypsometry_%s.png'%(self.figs_path,self.name_append))
    plt.savefig('%shypsometry_%s.png'%(self.figs_path,self.name_append))
