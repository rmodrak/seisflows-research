
from glob import glob
from os.path import basename, join

import numpy as np

from seisflows.tools import unix

from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system


class source_decimation_2d(loadclass('solver', 'source_decimation'), loadclass('solver', 'specfem2d')):
    """ Adds source_decimation machinery to SPECFEM2D
    """
    pass
