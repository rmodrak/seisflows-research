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


class GaussNewton(loadclass('postprocess', 'base')):

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

        # check input arguments
        if not exists(path):
            raise Exception()

        self.combine_kernels(path)
        self.process_kernels_lcg(path)

        # write gradient
        g = solver.merge(solver.load(
                PATH.HESS +'/'+ 'kernels/sum',
                verbose=True))

        g *= solver.merge(solver.load(
                path +'/'+ 'model'))

        # write gradient
        solver.save(PATH.HESS +'/'+ 'gauss_newton_lcg', solver.split(g))
        savenpy(PATH.OPTIMIZE +'/'+ 'g_lcg', g)


    def process_kernels_lcg(self, path):

        # process kernels
        if PAR.CLIP > 0.:
            system.run('solver', 'clip',
                       hosts='head',
                       path=path + '/' + 'kernels/sum',
                       thresh=PAR.CLIP)

        if PAR.SMOOTH > 0.:
            system.run('solver', 'smooth',
                       hosts='head',
                       path=path + '/' + 'kernels/sum',
                       span=PAR.SMOOTH)

