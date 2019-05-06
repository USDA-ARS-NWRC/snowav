
import matplotlib

def save(figure, paths):
    '''
    This function saves snowav figures.
    
    Args
        figure: figure object
        name: list of paths to save figure
    '''
    if type(paths) != list:
        paths = [paths]

    for path in paths:
        figure.savefig(path)
