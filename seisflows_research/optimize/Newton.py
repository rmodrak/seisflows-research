import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

from seisflows.optimize.lib.PLCG import PLCG

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


### utility functions


def _norm(v):
    return max(abs(v))



class Newton(custom_import('optimize', 'base')):
    """ Adds Newton-CG algorithm to nonlinear optimization base class
    """

    def check(cls):
        """ Checks parameters and paths
        """
        if 'SCHEME' not in  PAR:
            setattr(PAR, 'SCHEME', 'Newton')

        if 'LINESEARCH' not in PAR:
            setattr(PAR, 'LINESEARCH', 'Backtrack')

        if 'LCGPRECOND' not in PAR:
            setattr(PAR, 'LCGPRECOND', None)

        if 'LCGFORCE' not in PAR:
            setattr(PAR, 'LCGFORCE', 1.)

        if 'LCGMAX' not in PAR:
            setattr(PAR, 'LCGMAX', 2)

        if 'LCGTHRESH' not in PAR:
            setattr(PAR, 'LCGTHRESH', np.inf)

        if 'EPSILON' not in PAR:
            setattr(PAR, 'EPSILON', 1.)

        super(Newton, cls).check()


    def setup(cls):
        super(Newton, cls).setup()

        # prepare algorithm machinery
        cls.LCG = PLCG(
            path=PATH.OPTIMIZE, 
            thresh=PAR.LCGTHRESH, 
            maxiter=PAR.LCGMAX, 
            precond=PAR.LCGPRECOND,
            eta=PAR.LCGFORCE)


    def initialize_newton(cls):
        """ Initialize truncated Newton algorithm
        """
        unix.cd(PATH.OPTIMIZE)

        m = loadnpy('m_new')
        g = loadnpy('g_new')

        cls.LCG.initialize()

        cls.restarted = False

        v = loadnpy('LCG/p')
        h = PAR.EPSILON / _norm(v)

        savenpy('m_lcg', m + h*v)


    def iterate_newton(cls):
        """ Carry out truncated Newton iteration
        """
        unix.cd(PATH.OPTIMIZE)

        m = loadnpy('m_new')
        g = loadnpy('g_new')

        v = loadnpy('LCG/p')
        h = PAR.EPSILON / _norm(v)

        u = cls.hessian_product(g, h)

        isdone = cls.LCG.update(u)
        v = loadnpy('LCG/p')
        w = loadnpy('LCG/x')

        if not isdone:
            savenpy('m_lcg', m + h*v)
            return isdone

        s = np.dot(g,w)/np.dot(g,g)
        if s >= 0:
            print ' Newton direction rejected [not a descent direction]'
            w = -g
            s = 1.
            cls.restarted = True

        savenpy('p_new', w)
        savetxt('s_new', s)

        return isdone


    def hessian_product(cls, g, h):
        dg = loadnpy('g_lcg') - g
        return h**-1 * dg


    def initial_step(cls):
        return 1.


    def restart(cls):
        """ Discards history of algorithm; prepares to start again from 
          gradient direction
        """
        unix.cd(PATH.OPTIMIZE)

        g = cls.load('g_new')

        cls.save('p_new', -g)
        savetxt('s_new', cls.dot(g,g))

        if 'LBFGS' in PAR.LCGPRECOND:
            cls.LCG.LBFGS.restart()

        cls.restarted = True

        cls.stepwriter.iter -= 1
        cls.stepwriter.newline()


