
from copy import deepcopy
import numpy as np
import scipy.optimize

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.array import grid2mesh, mesh2grid, stack
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    loadclass, ParameterError

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import solver
import postprocess

migration = loadclass('workflow','migration')()


class compute_lambda(object):
    """ Postprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        if 'INPUT' not in PATH:
            setattr(PATH, 'INPUT', None)

        if not PATH.INPUT:
            migration.check()


    def main(self):
        """ Writes gradient of objective function
        """
        if not PATH.INPUT:
            migration.main()

        postprocess.process_kernels(
            path=PATH.GLOBAL)

        g = solver.load(PATH.GLOBAL +'/'+ 'kernels/sum', suffix='_kernel')
        m = solver.load(PATH.GLOBAL +'/'+ 'model')

        try:
            x = m['x'][0]
            z = m['z'][0]
            mesh = stack(x, z)
        except:
            from seisflows.seistools.io import loadbin
            model_path = PATH.OUTPUT +'/'+ 'model_true'
            x = loadbin(model_path, 0, 'x')
            z = loadbin(model_path, 0, 'z')
            mesh = stack(x, z)

        # apply mask
        if PATH.MASK:
            mask = solver.load(PATH.MASK)
            for key in solver.parameters:
                for iproc in range(PAR.NPROC):
                    g[key][iproc] *= mask[key][iproc]
        solver.save(PATH.OUTPUT +'/'+ 'gradient', g, suffix='_kernel')

        # compute spatial derivatives
        nabla = postprocess.nabla
        dg = {}
        for key in solver.parameters:
            dg[key] = []
            for iproc in range(PAR.NPROC):
                dg[key] += [nabla(mesh, m[key][iproc], g[key][iproc])]
        solver.save(PATH.OUTPUT +'/'+ 'nabla', dg, suffix='_kernel')

        # write regularized gradient
        kk = 0
        for lmbd in getlist(PAR.LAMBDA):
            gs = deepcopy(g)
            kk += 1
            for key in solver.parameters:
                for iproc in range(PAR.NPROC):
                    gs[key][iproc] += lmbd*dg[key][iproc]
            solver.save(PATH.OUTPUT +'/'+ 'gradient_%02d' % kk, gs, suffix='_kernel')

        # compute baseline
        lmbd0 = scipy.optimize.fmin_powell(
                 self.gradient_norm, 
                 0., 
                 ftol=1.e-9)
        print lmbd0


    def gradient_norm(self, lmbd, order=2, verbose=True):
        g = solver.merge(solver.load(
                PATH.OUTPUT +'/'+ 'gradient', suffix='_kernel'))

        dg = solver.merge(solver.load(
               PATH.OUTPUT +'/'+ 'nabla', suffix='_kernel'))

        f0 = np.linalg.norm(g, order)
        f = np.linalg.norm(g + lmbd*dg, order)

        if verbose:
            print ' %15.8e %15.8e' % (lmbd, f/f0)

        return f/f0


def getlist(var):
    if isinstance(var, list):
        return var

    if isinstance(var, float):
        return [var]

    if isinstance(var, int):
        return [var]
