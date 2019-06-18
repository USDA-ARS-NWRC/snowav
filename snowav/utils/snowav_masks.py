
import netCDF4 as nc
import os

def snowav_masks(dempath, plotorder = None, plotlabels = None, logger = None):
    '''
    Load dem, build snowav masks dictionary from topo.nc file. See
    CoreConfig.ini for more information on some config options and behavior.

    Args
    -----------
    dempath : str
        Path to topo.nc file, which is intended to be built by basin_setup
        package.
    plotorder : list
        Optional, default is all mask strings that appear in topo.nc file found
        in dempath.
    plotlabels : list
        Optional, default is mask_names.
    logger : list
        Optional, temporary snowav logger.

    Returns
    -----------
    out : dict
        dictionary keys are:
            dem: topo.nc dem
            nrows: number of rows in dem
            ncols: number of columns in dem
            masks: snowav required dictionary of masks found in topo.nc file
            plotorder: final version plotorder
            labels: dictionary with with plotorder[plotlabels] for optional
                different figure labels

    '''

    out = {}
    masks = {}
    labels = {}

    if logger is None:
        logger = []

    if not os.path.isfile(dempath):
        raise Exception('dempath {} not a valid file'.format(dempath))

    ncf = nc.Dataset(dempath, 'r')
    dem = ncf.variables['dem'][:]

    mask_names = [ncf.variables['mask'].long_name]

    for v in ncf.variables.items():

        # without the space before mask we get the basin total twice
        if ' mask' in v[0]:
            mask_names += [v[1].long_name]

    if plotorder is None:
        plotorder = mask_names

        # Exception for Cherry Creek
        plotorder =  ['Cherry Creek' if x=='Cherry' else x for x in plotorder]

    else:
        if len(plotorder) > len(mask_names):
            raise Exception(' Given [snowav] mask: {} '.format(plotorder) +
                            'to use with {}, but found '.format(dempath) +
                            '{}'.format(mask_names))

        for lbl in plotorder:
            if (lbl not in mask_names) and (lbl != 'Cherry Creek'):
                raise Exception(' Given [snowav] masks: {} '.format(lbl) +
                                'to use with {}, but found '.format(dempath) +
                                '{}'.format(mask_names))

    if plotlabels is None:
        for name in plotorder:
            labels[name] = name

    else:
        if len(plotlabels) != len(plotorder):
            logger.append(' Given [snowav] plotlabels ' +
                          '{} not equal to '.format(plotlabels) +
                          '{}, setting plotlabels '.format(plotorder) +
                          'to defaults')
            for name in plotorder:
                labels[name] = name

        else:
            for i, name in enumerate(plotorder):
                logger.append(' Assigning plot label {} '.format(plotlabels[i])+
                              'to mask {}'.format(name))
                labels[name] = plotlabels[i]

    nrows = len(dem[:,0])
    ncols = len(dem[0,:])

    for lbl in plotorder:

        # Exceptions
        if lbl == 'Cherry Creek':
            nclbl = 'Cherry'

        else:
            nclbl = lbl

        if lbl != plotorder[0]:
            masks[lbl] = {'mask':ncf[nclbl + ' mask'][:], 'label':lbl}

        else:
            masks[lbl] = {'mask':ncf['mask'][:], 'label':nclbl}

    ncf.close()

    out['dem'] = dem
    out['masks'] = masks
    out['nrows'] = nrows
    out['ncols'] = ncols
    out['plotorder'] = plotorder
    out['labels'] = labels
    out['logger'] = logger

    return out
