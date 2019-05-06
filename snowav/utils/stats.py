import numpy as np

def nashsutcliffe(observed,modeled):
    """
    Nash-Sutcliffe model efficinecy

        .. math::

         NSE = 1-\\frac{\\sum_{i=1}^{N}(e_{i}-s_{i})^2}{\\sum_{i=1}^{N}(e_{i}-\\bar{e})^2}

    :observed: Observed data
    :type: list

    :simulation: Modeled data
    :type: list

    :return: Nash-Sutcliff model efficiency
    :rtype: float

    """
    if len(observed) == len(modeled):
        m = np.array(modeled)
        o = np.array(observed)

        mean_observed = np.nanmean(o)

        numerator = np.nansum((m - o) ** 2)
        denominator = np.nansum((o - mean_observed)**2)

        return 1 - (numerator / denominator)

    else:
        logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan
