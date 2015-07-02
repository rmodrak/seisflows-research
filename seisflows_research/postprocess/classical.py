
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.array import grid2mesh, mesh2grid, nabla, stack
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver


class classical(loadclass('postprocess', 'base')):
    """ Adds classical regularization options to base class

        While the underlying theory is classical, strategies included here are
        also experimental in the sense that we are trying to optimize their
        performance on unstructured numerical grids.
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(Experimental, self).check()

        if 'REGULARIZE' not in PAR:
            setattr(PAR, 'REGULARIZE', None)

        if 'LAMBDA' not in PAR:
            setattr(PAR, 'LAMBDA', 0.)

        if 'MASK' not in PATH:
            raise ParameterError(PATH, 'MASK')

        if not exists(PATH.MASK):
            raise Exception


    def setup(self):
        """ Performs any required initialization or setup tasks
        """
        pass


    def process_kernels(self, path):
        """ Performs scaling, smoothing, and preconditioning operations in
            accordance with parameter settings
        """
        assert (exists(path))

        # apply smoothing
        if PAR.SMOOTH > 0.:
            system.run('solver', 'smooth',
                       hosts='head',
                       path=path + '/' + 'kernels/sum',
                       span=PAR.SMOOTH)

            # apply mask
            mask = solver.merge(solver.load(PATH.MASK))
            v_smooth = solver.merge(solver.load(path +'/'+ 'kernels/sum'))
            v_nosmooth = solver.merge(solver.load(path +'/'+ 'kernels/sum_nosmooth'))

            solver.save(path +'/'+ 'kernels/sum',
                        solver.split(v_smooth*mask),
                        suffix='_kernel')


    def regularize(self, path):
        assert (exists(path))

        # regularize
        if PAR.REGULARIZE:
            model = solver.load(path +'/'+ 'model')
            kernels = solver.load(path +'/'+ 'kernels/sum', suffix='_kernel')

            assert 'x' in model
            assert 'z' in model
            x = model['x'][0]
            z = model['z'][0]
            mesh = stack(x, z)

            for key in solver.parameters:            
                for iproc in range(PAR.NPROC):
                    kernel[key][iproc] += PAR.LAMBDA*
                        self.nabla, model[key][0], kernels[key][0])

            src = path +'/'+ 'kernels/sum'
            dst = path +'/'+ 'kernels/sum_noregularize'
            unix.mv(src, dst)

            solver.save(path +'/'+ 'kernels/sum',
                        kernels,
                        suffix='_kernel')


    def nabla(self, mesh, model, kernel):
        if PAR.REGULARIZE in ['TotalVariation']:
            V, grid = mesh2grid(kernel, mesh)
            DV = nabla(V, order=1)
            dv = grid2mesh(DV, grid, mesh)
            return np.sign(dv)/np.mean(model)

        elif PAR.REGULARIZE in ['Tikhonov0', 'Damping']:
            return v/np.mean(v)

        elif PAR.REGULARIZE in ['Tikhonov1', 'Smoothing']:
            V, grid = mesh2grid(v, mesh)
            DV = nabla(V, order=1)
            dv = grid2mesh(DV, grid, mesh)
            return dv/np.mean(v)

        elif PAR.REGULARIZE in ['Tikhonov2']:
            V, grid = mesh2grid(v, mesh)
            DV = nabla(V, order=2)
            dv = grid2mesh(DV, grid, mesh)
            return dv/np.mean(v)


