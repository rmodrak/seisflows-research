
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver


class Thomsen(loadclass('postprocess', 'base')):
    """ Postprocessing class
    """

    def write_gradient(self, path):
        """ Writes gradient of objective function
        """
        # check parameters
        if 'OPTIMIZE' not in PATH:
            raise ParameterError(PATH, 'OPTIMIZE')

        # check input arguments
        if not exists(path):
            raise Exception()

        self.combine_kernels(path)
        self.process_kernels(path)

        g = solver.merge(solver.load(
                path +'/'+ 'kernels/sum',
                suffix='_kernel',
                verbose=True))

        # apply scaling
        if float(PAR.SCALE) == 1.:
            pass
        elif not PAR.SCALE:
            pass
        else:
            g *= PAR.SCALE

        # write gradient
        solver.save(PATH.GRAD +'/'+ 'gradient', solver.split(g), suffix='_kernel')
        savenpy(PATH.OPTIMIZE +'/'+ 'g_new', g)


