
def plotlims(plotorder):
    '''
    Args
        basin: ['TUOL','SJ','KINGS','KAWEAH','MERCED','BRB','LAKES']
        plotorder: list of subbasin names (see database/tables.py)

    Returns
        lims: class with basin-specific plot details
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

        if basin == 'San Joaquin River Basin':
            pbbx = 0.3
            legy = 0.69
            legy2 = 0.65
            basin_name = 'sanjoaquin'
            # database_name = 'San Joaquin River Basin'

        if basin in ['RCEW']:
            pbbx = 0.1

        if basin == 'Kings River Basin':
            legy = 0.53
            legy2 = 0.39
            basin_name = 'kings'
            # database_name = 'Kings River Basin'

        if basin == 'Merced River Basin':
            legy = 0.58
            legy2 = 0.42
            basin_name = 'merced'
            # database_name = 'Merced River Basin'

        if basin == 'Kaweah River Basin':
            legy = 0.52
            legy2 = 0.38
            pbbx = 0.45
            basin_name = 'kaweah'
            # database_name = 'Kaweah River Basin'

        if basin == 'Lakes Basin':
            btx = 0.26
            legy = 0.85
            legy2 = 0.88
            basin_name = 'lakes'

        if basin == 'Extended Tuolumne':
            # legy = 0.74
            legy = 0.75
            btx = 0.28
            basin_name = 'tuolumne'

        if basin == 'Boise River Basin':
            # legy = 0.74
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
