
from os.path import join

from seisflows.seistools.io import copybin, loadbypar, savebin, splitvec, Minmax
from seisflows.seistools.io import Model as IOStruct

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()



class FwiAnisotropic2d(loadclass('solver', 'FwiElastic'), loadclass('solver', 'specfem2d_binary')):
    """ Adds elastic inversion machinery
    """
    model_parameters = []
    model_parameters += ['c11']
    model_parameters += ['c13']
    model_parameters += ['c15']
    model_parameters += ['c33']
    model_parameters += ['c35']
    model_parameters += ['c55']
    #model_parameters += ['c12']
    #model_parameters += ['c23']
    #model_parameters += ['c25']


    if PAR.MATERIALS == 'ChenTromp':
        from seisflows.seistools.maps import voigt_chentromp_2d as map_forward
        from seisflows.seistools.maps import chentromp_voigt_2d as map_inverse
        kernel_parameters = []
        kernel_parameters += ['A']
        kernel_parameters += ['C']
        kernel_parameters += ['N']
        kernel_parameters += ['L']
        kernel_parameters += ['F']


    elif PAR.MATERIALS == 'Voigt2d':
        from seisflows.seistools.maps import voigt_voigt_2d as map_forward
        from seisflows.seistools.maps import voigt_voigt_2d as map_inverse
        kernel_parameters = []
        kernel_parameters += ['c11']
        kernel_parameters += ['c13']
        kernel_parameters += ['c15']
        kernel_parameters += ['c33']
        kernel_parameters += ['c35']
        kernel_parameters += ['c55']


    elif PAR.MATERIALS == 'Thomsen2d':
        from seisflows.seistools.maps import voigt_thomsen_2d as map_forward
        from seisflows.seistools.maps import thomsen_voigt_2d as map_inverse
        kernel_parameters = []
        kernel_parameters += ['vp']
        kernel_parameters += ['vs']
        kernel_parameters += ['epsilon']
        kernel_parameters += ['delta']
        kernel_parameters += ['gamma']
        kernel_parameters += ['theta']

    else:
        raise ParameterError(PAR, 'Thomsen2d')



