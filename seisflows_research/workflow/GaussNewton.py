import numpy as np
import glob

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.tools import divides, exists
from seisflows.config import , \
    ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']
solver = sys.modules['seisflows_solver']
optimize = sys.modules['seisflows_optimize']
preprocess = sys.modules['seisflows_preprocess']
postprocess = sys.modules['seisflows_postprocess']


class GaussNewton(custom_import('workflow', 'Newton')):
    """ Inversion with truncated Gauss Newton model updates
    """

    def check_objects(self):
        if PAR.OPTIMIZE != 'GaussNewton':
            raise Exception

        if  PAR.PREPROCESS != 'GaussNewton':
            raise Exception

        if  PAR.POSTPROCESS != 'GaussNewton':
            raise Exception

