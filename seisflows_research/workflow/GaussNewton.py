import numpy as np
import glob

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import divides, exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver
import optimize
import preprocess
import postprocess


class GaussNewton(loadclass('workflow', 'Newton')):
    """ Inversion with truncated Gauss Newton model updates
    """

    def check_objects(self):
        if PAR.OPTIMIZE != 'GaussNewton':
            raise Exception

        if  PAR.PREPROCESS != 'GaussNewton':
            raise Exception

        if  PAR.POSTPROCESS != 'GaussNewton':
            raise Exception

