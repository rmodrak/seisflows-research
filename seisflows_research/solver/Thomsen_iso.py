
from seisflows.tools.config import loadclass


class Thomsen_iso(loadclass('solver', 'Thomsen_base')):

    # model parameters included in inversion
    parameters = []
    parameters += ['vp']
    parameters += ['vs']

