import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt
from seisflows.tools.config import loadclass, ParameterObj
from seisflows.optimize import lib

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class Newton(loadclass('optimize', 'base')):
    """ Implements truncated Newton algorithm
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
        p = cls.LCG.initialize()
        cls.delta = 1.e-3 * max(abs(p))**-1
        tmp = p*cls.delta
        savenpy('m_lcg', m + p*cls.delta)


    def iterate_newton(cls):
        """ Iterate truncated Newton algorithm
        """
        unix.cd(cls.path)

        g0 = loadnpy('g_new')
        g = loadnpy('g_lcg')

        p, isdone = cls.LCG.update((g - g0)/cls.delta)

        s = np.dot(g, p)
        if s >= 0:
            print ' Newton failed [not a descent direction]'
            p, isdone = -g, True

        if isdone:
            savenpy('p_new', p)
            savetxt('s_new', cls.delta*np.dot(g0, p))
        else:
            m = loadnpy('m_new')
            cls.delta = 1.e-3 * max(abs(p))**-1
            savenpy('m_lcg', m + p*cls.delta)

        return isdone

