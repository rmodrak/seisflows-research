
from os.path import join

from seisflows.seistools.io import copybin, loadbypar, savebin, splitvec, Minmax
from seisflows.seistools.io import Model as IOStruct

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()



class anisotropic2d(custom_import('solver', 'anisotropic'), custom_import('solver', 'specfem2d')):
    pass
