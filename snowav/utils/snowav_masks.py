

import netCDF4 as nc
import numpy as np
import os

def masks(dempath, plotorder = None, plotlabels = None):
    '''
    Load dem, build snowav masks dictionary from topo.nc file. See
    CoreConfig.ini for more information on some config options and behavior.

    Args
    ------
    dempath : str
        Path to topo.nc file, which is intended to be built by basin_setup
        package.
    plotorder : list
        Optional, default is all mask strings that appear in topo.nc file found
        in dempath.
    plotlabels : list
        Optional, default is mask_names.

    Returns
    ------
    out : dict

    '''

    out = {}
    masks = {}
    labels = {}
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
            raise Exception(' Config option [snowav] mask: {} for {}, but '
                            'found {}'.format(plotorder, dempath, mask_names))

        for lbl in plotorder:
            if (lbl not in mask_names) and (lbl != 'Cherry Creek'):
                raise Exception(' Config option [snowav] mask: {} for {}, but '
                                'found {}'.format(plotorder, dempath, mask_names))

    if plotlabels is None:
        for name in plotorder:
            labels[name] = name

    else:
        if len(plotlabels) != len(plotorder):
            logger.append(' Config option [snowav] plotlabels {} not equal to '
                          '{}, setting plotlabels to '
                          'defaults'.format(plotlabels, plotorder))
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

def precip(hr, rundirs_dict, path):
    '''
    Get daily total precip and percent rain from hourly smrf outputs in
    expected awsm_daily format.

    Args
    -----
    hr : int
    rundirs_dict : dict
    path : str

    Returns
    -------
    flag : bool
    path : str
    rain : array
    precip : array

    '''

    sf = rundirs_dict[hr].replace('runs','data')
    sf = sf.replace('run','data')
    ppt_path = os.path.joint(sf,path)
    percent_snow_path = ppt_path.replace('precip','percent_snow')

    # precip = np.zeros((args['nrows'],args['ncols']))
    # rain = np.zeros((args['nrows'],args['ncols']))

    if os.path.isfile(ppt_path):
        ppt = nc.Dataset(ppt_path, 'r')
        percent_snow = nc.Dataset(percent_snow_path, 'r')

        # For the wy2019 daily runs, precip.nc always has an extra hour, but
        # in some WRF forecast runs there are fewer than 24...
        if len(ppt.variables['precip'][:]) > 24:
            nb = 24

        else:
            nb = len(ppt.variables['precip'][:])

        for nb in range(0,nb):
            pre = ppt['precip'][nb]
            ps = percent_snow['percent_snow'][nb]

            if nb == 0:
                rain = np.multiply(pre,(1-ps))
                precip = copy.deepcopy(pre)

            else:
                rain = rain + np.multiply(pre,(1-ps))
                precip = precip + copy.deepcopy(pre)

        ppt.close()
        percent_snow.close()

    return flag, ppt_path, rain, precip
