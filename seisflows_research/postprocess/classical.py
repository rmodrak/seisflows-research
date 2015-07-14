
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.array import grid2mesh, mesh2grid, stack
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass
from seisflows.tools.math import nabla


PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver


class classical(loadclass('postprocess', 'base')):
    """ Adds classical regularization options to base class

        While the underlying theory is completely classical, some options are
        experimental in the sense that we are trying to improve performance on
        unstructured numerical grids.
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(classical, self).check()

        if 'RADIUS' not in PAR:
            raise ParameterError(PAR, 'RADIUS')

        if 'REGULARIZE' not in PAR:
            setattr(PAR, 'REGULARIZE', None)

        if 'LAMBDA' not in PAR:
            setattr(PAR, 'LAMBDA', 0.)

        # assertions
        assert PAR.REGULARIZE in ['', None,
                'TotalVariation',
                'Tikhonov0', 'Damping',
                'Tikhonov1', 'Smoothing',
                'Tikhonov2']

    def setup(self):
        """ Performs any required initialization or setup tasks
        """
        pass


    def process_kernels(self, path, debug=False):
        """ Masks source and receiver artifacts
        """
        fullpath = path +'/'+ 'kernels'
        assert exists(path)

        if debug:
            unix.cp(fullpath, fullpath+'_nomask')

            system.run('solver', 'combine',
                       hosts='head',
                       path=fullpath+'_nomask')

        # mask sources and receivers
        system.run('postprocess', 'mask', 
                   hosts='all', 
                   path=fullpath)

        system.run('solver', 'combine',
                   hosts='head',
                   path=fullpath)


    def test_regularize(self):
        from copy import deepcopy

        path = PATH.INPUT

        kernels = solver.load(path +'/'+ 'kernels/sum')
        if PATH.MASK:
            mask = solver.load(PATH.MASK)
            for key in solver.parameters:
                for iproc in range(PAR.NPROC):
                    kernels[key][iproc] *= mask[key][iproc]


        model = solver.load(path +'/'+ 'model')
        assert 'x' in model
        assert 'z' in model
        x = model['x'][0]
        z = model['z'][0]
        mesh = stack(x, z)

        for Lambda in np.linspace(0., PAR.NN*PAR.HH, PAR.NN+1):
            output = deepcopy(kernels)
            print output.keys()
            for key in solver.parameters:
                for iproc in range(PAR.NPROC):
                    output[key][iproc] += Lambda * \
                        self.nabla(mesh, model[key][iproc], kernels[key][iproc])

            solver.save(PATH.OUTPUT +'/'+ 'gradient_%04d' % int(1.e-2*Lambda),
                        output,
                        suffix='_kernel')


    def regularize(self, path):
        assert (exists(path))

        # regularize
        if PAR.REGULARIZE:
            model = solver.load(path +'/'+ 'model')
            kernels = solver.load(path +'/'+ 'gradient', suffix='_kernel')

            assert 'x' in model
            assert 'z' in model
            x = model['x'][0]
            z = model['z'][0]
            mesh = stack(x, z)

            for key in solver.parameters:            
                for iproc in range(PAR.NPROC):
                    kernels[key][iproc] +=\
                        PAR.LAMBDA * self.nabla(mesh, model[key][iproc], kernels[key][iproc])

            return solver.merge(kernels)


    def nabla(self, mesh, model, kernel):
        if PAR.REGULARIZE in ['TotalVariation']:
            v = model
            V, grid = mesh2grid(kernel, mesh)
            DV = nabla(V, order=1)
            dv = grid2mesh(DV, grid, mesh)
            return -np.sign(dv)/np.mean(v)

        elif PAR.REGULARIZE in ['Tikhonov0', 'Damping']:
            v = model
            return v/np.mean(v)

        elif PAR.REGULARIZE in ['Tikhonov1', 'Smoothing']:
            v = model
            V, grid = mesh2grid(kernel, mesh)
            DV = nabla(V, order=1)
            dv = grid2mesh(DV, grid, mesh)
            return dv/np.mean(v)

        elif PAR.REGULARIZE in ['Tikhonov2']:
            v = model
            V, grid = mesh2grid(kernel, mesh)
            DV = nabla(V, order=2)
            dv = grid2mesh(DV, grid, mesh)
            return dv/np.mean(v)


    def mask(self, path=''):
        import preprocess

        fullpath = path +'/'+  solver.getname
        kernels = solver.load(fullpath)
        if not PAR.RADIUS:
            return

        x = kernels['x'][0]
        z = kernels['z'][0]
        mesh = stack(x, z)
        _, h = preprocess.load(solver.getpath +'/'+ 'traces/obs')

        # mask sources
        mask = np.exp(-0.5*((x-h.sx[0])**2.+(z-h.sy[0])**2.)/PAR.RADIUS**2.)
        for key in solver.parameters:
            weight = np.sum(mask*kernels[key][0])/np.sum(mask)
            kernels[key][0] *= 1.-mask
            kernels[key][0] += mask*weight

        # mask receivers
        for ir in range(h.nr):
            mask = np.exp(-0.5*((x-h.rx[ir])**2.+(z-h.ry[ir])**2.)/PAR.RADIUS**2.)
            for key in solver.parameters:
                weight = np.sum(mask*kernels[key][0])/np.sum(mask)
                kernels[key][0] *= 1.-mask
                kernels[key][0] += mask*weight

        solver.save(fullpath, kernels)

