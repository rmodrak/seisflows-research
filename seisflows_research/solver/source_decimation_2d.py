
from glob import glob
from os.path import basename, join

import numpy as np

from seisflows.tools import unix

from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system


class source_decimation_2d(custom_import('solver', 'source_decimation'), custom_import('solver', 'specfem2d')):
    """ Adds source_decimation machinery to SPECFEM2D
    """
    pass
