
from glob import glob
from os.path import join

from seisflows.tools import unix

from seisflows.seistools.io import copybin, savebin
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system


class Thomsen_base(loadclass('solver', 'specfem3d_legacy')):

    raise NotImplementedError("The following methods need to be fixed: export_kernels")

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
            parameters = self.solver_parameters
            parameters += ['rho']
            for parameter in parameters:
                unix.mkdir(path)
                src = glob(join(self.model_databases, '*'+parameter+'.bin'))
                dst = path
                unix.cp(src, dst)


    def export_kernels(self, path):
        raise NotImplementedError
        try:
            files = glob(self.model_databases +'/'+ '*alpha*_kernel.bin')
            unix.rename('alpha', 'vp', files)

            files = glob(self.model_databases +'/'+ '*beta*_kernel.bin')
            unix.rename('beta', 'vs', files)
        except:
            pass

        # export kernels
        unix.mkdir_gpfs(join(path, 'kernels'))
        unix.mkdir(join(path, 'kernels', self.getname))
        src = join(glob(self.model_databases +'/'+ '*kernel.bin'))
        dst = join(path, 'kernels', self.getname)
        unix.mv(src, dst)


