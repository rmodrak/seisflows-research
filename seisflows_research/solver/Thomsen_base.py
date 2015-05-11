
from glob import glob
from os.path import join

from seisflows.tools import unix

from seisflows.seistools.io import copybin, savebin
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system


class Thomsen_base(loadclass('solver', 'specfem3d_legacy')):

    #raise NotImplementedError("The following methods need to be fixed: export_kernels")

    # parameters expected by solver
    solver_parameters = []
    solver_parameters += ['vp']
    solver_parameters += ['vs']
    solver_parameters += ['epsilon']
    solver_parameters += ['delta']
    solver_parameters += ['gamma']
    solver_parameters += ['theta']
    solver_parameters += ['azimuth']


    def save(self, path, model, prefix='', suffix=''):
        super(Thomsen_base, self).save(
            path, model, prefix, suffix, self.solver_parameters)


    def export_model(self, path):
        super(Thomsen_base, self).export_model(
            path, self.solver_parameters+['rho'])


