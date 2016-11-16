
from os.path import join

import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.array import grid2mesh, mesh2grid, stack
from seisflows.tools.code import exists, Struct
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import
from seisflows.tools.math import gauss2

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver


class model_projection(custom_import('postprocess', 'base')):
    """ Projects from GLL mesh to Gaussian function basis

        SO FAR, CAN ONLY BE USED FOR ACOUSTIC 2D WAVEFORM INVERSION.
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(model_projection, self).check()

        if 'REFERENCE' not in PAR:
            setattr(PAR, 'REFERENCE', None)

        assert PAR.SOLVER in ['specfem2d', 'specfem2d_legacy']
        assert PAR.MATERIALS in ['Acoustic']


    def setup(self):
        """ Sets up basis functions
        """
        super(model_projection, self).setup()

        x,z = self.getxz()
        dx = (x.max()-x.min())/(PAR.NX+2)
        dz = (z.max()-z.min())/(PAR.NZ+2)

        bx = np.linspace(x.min()+dx/2, x.max()-dx/2, (PAR.NX))
        bz = np.linspace(z.min()+dz/2, z.max()-dz/2, (PAR.NZ))

        self.nodes = Struct()
        self.nodes.x = x
        self.nodes.z = z
        self.nodes.dx = dx
        self.nodes.dz = dz
        self.nodes.bx = bx
        self.nodes.bz = bz

        c = self.project_from_gll(np.ones(x.shape))
        w = self.project_to_gll(c, weight=False)**-1
        self.nodes.weights = w

        if exists(PATH.MODEL_TRUE):
            # write true model
            m = solver.merge(solver.load(PATH.MODEL_TRUE))
            c = self.project_from_gll(m)
            m = self.project_to_gll(c)
            path = join(PATH.OUTPUT, 'model_true_projected')
            solver.save(path, solver.split(m))

        # write initial model
        m = solver.merge(solver.load(PATH.MODEL_INIT))
        c = self.project_from_gll(m)
        m = self.project_to_gll(c)
        path = join(PATH.OUTPUT, 'model_init_projected')
        solver.save(path, solver.split(m))

        # overwrite m_new
        unix.mkdir(PATH.OPTIMIZE)
        if PAR.REFERENCE:
            filename = join(PATH.OPTIMIZE, 'm_new')
            savenpy(filename, np.zeros(PAR.NX*PAR.NZ))
        else:
            filename = join(PATH.OPTIMIZE, 'm_new')
            savenpy(filename, c)


    def write_gradient(self, path):
        """ Reads kernels and writes gradient of objective function
        """
        if 'OPTIMIZE' not in PATH:
            raise ParameterError(PATH, 'OPTIMIZE')

        if not exists(path):
            raise Exception()

        self.combine_kernels(
            path=path,
            parameters=solver.parameters)

        g = solver.merge(solver.load(
                 path +'/'+ 'kernels/sum',
                 suffix='_kernel',
                 verbose=True))

        if PAR.LOGARITHMIC:
            # convert from logarithmic to absolute perturbations
            g *= solver.merge(solver.load(path +'/'+ 'model'))
        self.save(path, g)

        if PATH.MASK:
            # apply mask
            g *= solver.merge(solver.load(PATH.MASK))
            self.save(path, g, backup='nomask')

        c = self.project_from_gll(g)
        savenpy(PATH.OPTIMIZE +'/'+ 'g_new', c)


    def project_from_gll(self, v):
        x = self.nodes.x
        z = self.nodes.z
        c = []
        for ib in range((PAR.NX*PAR.NZ)):
           c += [np.dot(self.gauss2(ib,x,z), v)]
        return np.array(c)


    def project_to_gll(self, c, weight=True):
        x = self.nodes.x
        z = self.nodes.z
        v = np.zeros(x.shape)

        for ib in range((PAR.NX*PAR.NZ)):
            v += c[ib] * self.gauss2(ib,x,z)[:]

        if weight:
            v *= self.nodes.weights
        return v


    def getxz(self):
        model_path = PATH.MODEL_TRUE
        try:
            m = solver.load(model_path)
            x = m['x'][0]
            z = m['z'][0]
        except:
            from seisflows.plugins.io import loadbin
            x = loadbin(model_path, 0, 'x')
            z = loadbin(model_path, 0, 'z')
        return x,z


    def gauss2(self, ib,x,z):
        ix = ib % PAR.NX
        iz = ib / PAR.NX

        dx = self.nodes.dx
        dz = self.nodes.dz
        bx = self.nodes.bx
        bz = self.nodes.bz

        mu = [bx[ix], bz[iz]]
        sigma = np.array([[dx**2,0], [0,dz**2]])

        return gauss2(x,z,mu,sigma, normalize=False)


#   def process_kernels(self, path, parameters):
#        """ Processes kernels in accordance with parameter settings
#        """
#        fullpath = path +'/'+ 'kernels'
#        assert exists(path)
#
#        if exists(fullpath +'/'+ 'sum'):
#            unix.mv(fullpath +'/'+ 'sum', fullpath +'/'+ 'sum_nofix')
#
#        # mask sources and receivers
#        system.run('postprocess', 'fix_near_field',
#                   hosts='all',
#                   path=fullpath)
#
#        system.run('solver', 'combine',
#                   hosts='head',
#                   path=fullpath,
#                   parameters=parameters)

