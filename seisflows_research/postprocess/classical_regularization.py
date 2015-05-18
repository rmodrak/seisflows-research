
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


class classical_regularization(loadclass('postprocess', 'base')):
    """ Adds alternate regularization options to base class

        Regularization strategies are experimental in the sense that 
        we are trying to make them perform better when the model is 
        expressed on an unstructured numerical grid.
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

        if PATH.MASK:
            assert exists(PATH.MASK)


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

        try:
            v = solver.merge(solver.load(path +'/'+ 'kernels/sum_nosmooth'))
            w = solver.merge(solver.load(path +'/'+ 'kernels/sum'))
        except:
            v = solver.merge(solver.load(path +'/'+ 'kernels/sum'))
            w = v

        solver.save(path +'/'+ 'kernels/sum',
                    solver.split(v*(mask) + w*(1.-mask)),
                    suffix='_kernel')


    def regularize(self, path):
        """ Performs scaling, smoothing, and preconditioning operations in
            accordance with parameter settings
        """
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


