
def plotlims(plotorder):
    '''
    Compile basin-specific limits for making tidy figures.

    Args
    ------
    plotorder : list
        subbasin names

    Returns
    ------
    lims : class
        basin-specific plot details:
        pbbx:
        legx: legend x
        legy: legend y
        btx: basin total text x
        bty: basin total text y
        sumorder: list of subbasin to plot (omits total unless single basin)
        swid: bar spacing

    '''
    basin = plotorder[0]

    class lims:

        # snow-free patch label
        pbbx = 0.05

        # legend x,y
        legx = 0.01
        legy = 0.71

        # second legend x,y for precip_depth
        legy2 = 0.71

        btx = 0.31
        bty = 0.95

        if 'San Joaquin' in basin:
            pbbx = 0.3
            legy = 0.69
            legy2 = 0.65
            basin_name = 'sanjoaquin'

        if basin in ['RCEW']:
            pbbx = 0.1

        if 'Kings' in basin:
            legy = 0.52
            legy2 = 0.39
            basin_name = 'kings'

        if 'Merced' in basin:
            legy = 0.58
            legy2 = 0.42
            basin_name = 'merced'

        if 'Kaweah' in basin:
            legy = 0.52
            legy2 = 0.38
            pbbx = 0.45
            basin_name = 'kaweah'

        if 'Lakes' in basin:
            btx = 0.26
            legy = 0.85
            legy2 = 0.88
            basin_name = 'lakes'

        if basin == 'Extended Tuolumne':
            legy = 0.75
            btx = 0.28
            basin_name = 'tuolumne'

        if basin == 'Tuolumne River Basin':
            legy = 0.46
            btx = 0.28
            pbbx = 0.25
            legy2 = 0.30
            basin_name = 'tuolumne'

        if 'Boise' in basin:
            legy = 0.75
            btx = 0.28
            basin_name = 'brb'

        # Total basin label
        if plotorder is not None:
            if len(plotorder) > 1:
                sumorder = plotorder[1::]
                swid = 0.25

            else:
                sumorder = plotorder
                swid = 0.45

    return lims
