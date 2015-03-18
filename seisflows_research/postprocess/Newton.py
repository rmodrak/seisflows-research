import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists
from seisflows.tools.config import loadclass, ParameterObj
from seisflows.optimize import lib

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import solver
import system


class Newton(loadclass('postprocess', 'base')):

    def check(self):
        super(Newton, self).check()

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
                type='model',
                verbose=True))

        g *= solver.merge(solver.load(
                path +'/'+ 'model'))

        savenpy(PATH.OPTIMIZE +'/'+ 'g_lcg', g)

        # write gradient
        solver.save(PATH.HESS +'/'+ 'newton_lcg', solver.split(g))
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

