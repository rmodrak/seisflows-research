import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

from seisflows.optimize import lib

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class GaussNewton(loadclass('optimize', 'Newton')):
    """ Adds Gauss-Newton-CG algorithm to nonlinear optimization base class
    """

    def check(cls):
        """ Checks parameters and paths
        """
        super(GaussNewton, cls).check()

        if 'SCHEME' not in  PAR:
            setattr(PAR, 'SCHEME', 'GaussNewton')

        if 'LCGMAX' not in PAR:
            setattr(PAR, 'LCGMAX', 2)

        if 'LCGTHRESH' not in PAR:
            setattr(PAR, 'LCGTHRESH', np.inf)

        if 'LCGPRECOND' not in PAR:
            setattr(PAR, 'LCGPRECOND', 0)


    def hessian_product(cls, g, h):
        unix.cd(cls.path)

        dg = loadnpy('g_lcg')
        return h**-1 * dg


