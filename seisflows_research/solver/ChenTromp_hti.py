from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class ChenTromp_hti(loadclass('solver', 'ChenTromp_base')):

    # model parameters included in inversion
    parameters = []
    parameters += ['A']
    parameters += ['C']
    parameters += ['L']
    parameters += ['N']
    parameters += ['F']
    parameters += ['Gc']
    parameters += ['Gs']
    parameters += ['Bc']
    parameters += ['Bs']
    parameters += ['Hc']
    parameters += ['Hs']
    parameters += ['Ec']
    parameters += ['Es']

