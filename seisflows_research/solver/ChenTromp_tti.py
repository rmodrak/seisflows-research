from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class ChenTromp_tti(loadclass('solver', 'ChenTromp_base')):

    # model parameters included in inversion
    parameters = []
    parameters += ['A']
    parameters += ['C']
    parameters += ['L']
    parameters += ['N']
    parameters += ['F']
    parameters += ['Jc']
    parameters += ['Js']
    parameters += ['Kc']
    parameters += ['Ks']
    parameters += ['Mc']
    parameters += ['Ms']
    parameters += ['Gc']
    parameters += ['Gs']
    parameters += ['Bc']
    parameters += ['Bs']
    parameters += ['Hc']
    parameters += ['Hs']
    parameters += ['Dc']
    parameters += ['Ds']
    parameters += ['Ec']
    parameters += ['Es']

