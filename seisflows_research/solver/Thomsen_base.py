
from glob import glob
from os.path import join

from seisflows.tools import unix

from seisflows.seistools.io import copybin, savebin
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system


class Thomsen_base(loadclass('solver', 'specfem3d_legacy')):
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
        unix.mkdir(path)

        for iproc in range(PAR.NPROC):
            for key in self.solver_parameters:
                if key in self.parameters:
                    savebin(model[key][iproc], path, iproc, prefix+key+suffix)
                elif 'kernel' not in suffix:
                    src = PATH.OUTPUT +'/'+ 'model_init'
                    dst = path
                    copybin(src, dst, iproc, prefix+key+suffix)

            if 'rho' in self.parameters:
                savebin(model['rho'][iproc], path, iproc, prefix+'rho'+suffix)
            else:
                src = PATH.OUTPUT +'/'+ 'model_init'
                dst = path
                copybin(src, dst, iproc, 'rho')


    def export_model(self, path):
        if system.getnode() == 0:
            for parameter in self.solver_parameters:
                unix.mkdir(path)
                src = glob(join(self.model_databases, '*'+parameter+'.bin'))
                dst = path
                unix.cp(src, dst)

