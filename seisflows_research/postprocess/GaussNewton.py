import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import solver
import system


class base(loadclass('postprocess', 'base')):

    def check(self):
        super(GaussNewton, self).check()

        # check paths
        if 'OPTIMIZE' not in PATH:
            raise ParameterError(PATH, 'OPTIMIZE')

        if 'GRAD' not in PATH:
            raise ParameterError(PATH, 'GRAD')

        if 'HESS' not in PATH:
            raise ParameterError(PATH, 'HESS')


    def write_gradient_lcg(self, path):
        """ Reads kernels and writes gradient of objective function
        """
        if 'OPTIMIZE' not in PATH:
            raise ParameterError(PATH, 'OPTIMIZE')

        if not exists(path):
            raise Exception()

        self.combine_kernels(path, solver.parameters)
        self.process_kernels(path, solver.parameters)

        g = solver.merge(solver.load(
                 path +'/'+ 'kernels/sum',
                 suffix='_kernel',
                 verbose=True))

        if PAR.LOGARITHMIC:
            # convert from logarithmic to absolute perturbations
            g *= solver.merge(solver.load(path +'/'+ 'model'))
        self.save(path, g)

        if PATH.MASK:
            # apply mask
            g *= solver.merge(solver.load(PATH.MASK))
            self.save(path, g, backup='nomask')

        # write gradient
        solver.save(PATH.HESS +'/'+ 'gauss_newton_lcg', solver.split(g))
        savenpy(PATH.OPTIMIZE +'/'+ 'g_lcg', g)


GaussNewton = base

