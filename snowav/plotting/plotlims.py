
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

        # snow-free patch label
        pbbx = 0.05

        # legend x,y
        legx = 0.01
        legy = 0.71

        # second legend x,y for precip_depth
        legy2 = 0.71

        btx = 0.31
        bty = 0.95

        if basin in ['SJ']:
            pbbx = 0.3
            legy = 0.69
            legy2 = 0.65

        if basin in ['RCEW']:
            pbbx = 0.1

        if basin in ['KINGS']:
            legy = 0.52
            legy2 = 0.39

        if basin in ['MERCED']:
            legy = 0.57
            legy2 = 0.42

        if basin in ['KAWEAH']:
            legy = 0.5
            legy2 = 0.39
            pbbx = 0.45

        if basin in ['LAKES']:
            btx = 0.26
            legy = 0.85
            legy2 = 0.88

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

        # if ['Main'] in sumorder:
        #     sumorder.replace('Main','Mammoth')

    return lims
