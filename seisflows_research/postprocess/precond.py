
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


class precond(loadclass('postprocess', 'base')):
    """ Adds additional preconditioners to base class
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(debug, self).check()        


    def setup(self):
        """ Performs any required initialization or setup tasks
        """
        super(debug, self).setup()

        if PAR.PRECOND.lower() in ['pca', 'principal_components']:
            self.precond = principle_components

        elif PAR.PRECOND.lower() in ['diag', 'diagonal', 'user_supplied_diagonal']:
            assert exists(PATH.PRECOND)
            self.precond = user_supplied_diagonal


### preconditioners

def user_supplied_diagonal(q):
    p = solver.merge(solver.load(PATH.PRECOND))
    return p*q

def principal_components(q):
    raise NotImplementedError
    r = solver.split(q)

    # compute covariance
    c = {}
    for key1 in solver.parameters:
        c[key1] = {}
        for key2 in solver.parameters:
            c[key1][key2] = np.dot(r[key1], r[key2])

    # diagonalize
    a = eig(c)
    b = inv(a)

    # apply preconditioner
    s = {}
    for key1 in solver.parameters:
        for key2 in solver.parameters:
            for iproc in range(solver.mesh.nproc):
                s[key1][iproc] +=  [d[key2]*r[key][iproc]]

