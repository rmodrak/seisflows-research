
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
    model_parameters += ['C11']
    model_parameters += ['C13']
    model_parameters += ['C15']
    model_parameters += ['C33']
    model_parameters += ['C35']
    model_parameters += ['C55']
    model_parameters += ['C12']
    model_parameters += ['C23']
    model_parameters += ['C25']


    if PAR.MATERIALS == 'ChenTromp':
        from seisflows.seistools.maps import ChenTromp2d_forward as map_forward
        from seisflows.seistools.maps import ChenTromp2d_forward as map_inverse
        kernel_parameters = []
        kernel_parameters += ['A']
        kernel_parameters += ['C']
        kernel_parameters += ['N']
        kernel_parameters += ['L']
        kernel_parameters += ['F']


    elif PAR.MATERIALS == 'Voigt2d':
        from seisflows.seistools.maps import Voigt2d_forward as map_forward
        from seisflows.seistools.maps import Voigt2d_inverse as map_inverse
        kernel_parameters = []
        kernel_parameters += ['C11']
        kernel_parameters += ['C13']
        kernel_parameters += ['C15']
        kernel_parameters += ['C33']
        kernel_parameters += ['C35']
        kernel_parameters += ['C55']
        kernel_parameters += ['C12']
        kernel_parameters += ['C23']
        kernel_parameters += ['C25']


    elif PAR.MATERIALS == 'Thomsen2d':
        from seisflows.seistools.maps import Thomsen2d_forward as map_forward
        from seisflows.seistools.maps import Thomsen2d_inverse as map_inverse
        kernel_parameters = []
        kernel_parameters += ['vp']
        kernel_parameters += ['vs']
        kernel_parameters += ['epsilon']
        kernel_parameters += ['delta']
        kernel_parameters += ['gamma']
        kernel_parameters += ['theta']



