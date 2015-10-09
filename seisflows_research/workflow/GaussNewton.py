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


class GaussNewton(loadclass('workflow', 'inversion')):
    """ Inversion with truncated Gauss Newton model updates
    """

    def check(self):
        super(GaussNewton, self).check()

        if PAR.OPTIMIZE != 'GaussNewton':
            raise Exception

        if  PAR.PREPROCESS != 'GaussNewton':
            raise Exception

        if  PAR.POSTPROCESS != 'GaussNewton':
            raise Exception


    def compute_direction(self):
        """ Computes search direction
        """
        self.evaluate_gradient()

        optimize.initialize_newton()
        for self.ilcg in range(1, PAR.LCGMAX+1):
            self.apply_hessian()
            isdone = optimize.iterate_newton()
            if isdone:
                break


    def apply_hessian(self):
        """ Computes the action of the Hessian on a given vector through
          solver calls
        """
        if PAR.VERBOSE:
            print " LCG iteration", self.ilcg

        self.write_model(path=PATH.HESS, suffix='lcg')

        system.run('solver', 'apply_hess',
                   hosts='all',
                   path=PATH.HESS)

        postprocess.write_gradient_lcg(
            path=PATH.HESS)

        unix.rm(PATH.HESS+'_debug')
        unix.mv(PATH.HESS, PATH.HESS+'_debug')
        unix.mkdir(PATH.HESS)



