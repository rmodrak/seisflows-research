import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

from seisflows.optimize import lib

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


### utility functions

_eps = 1.e-3

def _norm(v):
    return max(abs(v))



class Newton(loadclass('optimize', 'base')):
    """ Adds Newton-CG algorithm to nonlinear optimization base class
    """

    def check(cls):
        """ Checks parameters and paths
        """
        super(Newton, cls).check()

        if 'SCHEME' not in  PAR:
            setattr(PAR, 'SCHEME', 'Newton')

        if 'LCGMAX' not in PAR:
            setattr(PAR, 'LCGMAX', 2)

        if 'LCGTHRESH' not in PAR:
            setattr(PAR, 'LCGTHRESH', np.inf)

        if 'LCGPRECOND' not in PAR:
            setattr(PAR, 'LCGPRECOND', 0)


    def setup(cls):
        super(Newton, cls).setup()
        cls.LCG = lib.LCG(cls.path, PAR.LCGTHRESH, PAR.LCGMAX, PAR.LCGPRECOND)


    def initialize_newton(cls):
        """ Initialize truncated Newton algorithm
        """
        unix.cd(cls.path)

        m = loadnpy('m_new')
        g = loadnpy('g_new')

        h = _eps / _norm(g)

        cls.LCG.initialize()
        v = loadnpy('LCG/p')

        savenpy('m_lcg', m + h*v)


    def iterate_newton(cls):
        """ Carry out truncated Newton iteration
        """
        unix.cd(cls.path)

        m = loadnpy('m_new')
        g = loadnpy('g_new')

        h = _eps / _norm(g)
        u = cls.hessian_product(g, h)

        isdone = cls.LCG.update(u)
        v = loadnpy('LCG/p')
        w = loadnpy('LCG/x')
        s = np.dot(g,w)/np.dot(g,g)

        if not isdone:
            savenpy('m_lcg', m + h*v)
            return isdone

        if s >= 0:
            print ' Stopping Newton algorithm [not a descent direction]'
            w = -g
            s = 1.

        savenpy('p_new', w)
        savetxt('s_new', s)

        return isdone


    def hessian_product(cls, g, h):
        dg = loadnpy('g_lcg') - g
        return h**-1 * dg


