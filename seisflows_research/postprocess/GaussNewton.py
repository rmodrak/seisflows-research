import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt
from seisflows.tools.config import loadclass, ParameterObj
from seisflows.optimize import lib

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import solver


class GaussNewton(loadclass('postprocess', 'base')):

    def process_kernels(self, *args, **kwargs):
        super(GaussNewton, self).process_kernels(*args, **kwargs)

        if 'tag' not in kwargs:
            return

        if kwargs['tag'] != 'gauss_newton_lcg':
            return

        g = solver.merge(solver.load(
                PATH.HESS +'/'+ 'gauss_newton_lcg', 
                type='model', 
                verbose=True))

        g0 = solver.merge(solver.load(
                 PATH.GRAD +'/'+ 'gradient',
                 type='model', 
                 verbose=True))

        savenpy(PATH.OPTIMIZE +'/'+ 'g_lcg', g)


        if 0:
           # print debugging info
           dg = g-g0
           print ' %e  %e' % (g.min(), g.max())
           print ' %e  %e' % (g0.min(), g0.max())
           print ' %e  %e' % (dg.min(), dg.max())

