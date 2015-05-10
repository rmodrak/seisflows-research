
from seisflows.tools.config import loadclass


class ChenTromp_vti(loadclass('solver', 'ChenTromp_base')):

    # model parameters included in inversion
    parameters = []
    parameters += ['A']
    parameters += ['C']
    parameters += ['L']
    parameters += ['N']
    parameters += ['F']

