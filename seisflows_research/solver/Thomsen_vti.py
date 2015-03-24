from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class Thomsen_vti(loadclass('solver', 'Thomsen_base')):

    # model parameters included in inversion
    inversion_parameters = []
    inversion_parameters += ['vp']
    inversion_parameters += ['vs']
    inversion_parameters += ['epsilon']
    inversion_parameters += ['delta']
    inversion_parameters += ['gamma']

