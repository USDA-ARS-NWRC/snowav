
def plotlims(basin,plotorder):
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

    class lims:

        pbbx = 0.05
        legx = 0.01
        legy = 0.71
        btx = 0.31
        bty = 0.95

        if basin == 'SJ':
            pbbx = 0.3
            legy = 0.69

        if basin in [ 'RCEW']:
            pbbx = 0.1

        if basin in ['KINGS', 'MERCED', 'KAWEAH']:
            legy = 0.5

        if basin in ['LAKES']:
            btx = 0.26

        if basin in ['BRB','TUOL']:
            legy = 0.74
            btx = 0.28

        # Total basin label
        if len(plotorder) > 1:
            sumorder = plotorder[1::]
            swid = 0.25
        else:
            sumorder = plotorder
            swid = 0.45

    return lims
