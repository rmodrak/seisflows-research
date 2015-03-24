
from seisflows.tools import unix

from seisflows.seistools.io import copybin, savebin
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class Thomsen_base(loadclass('solver', 'specfem3d_legacy')):
    # parameters expected by solver
    solver_parameters = []
    solver_parameters += ['rho']
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
            for key in solver.parameters:
                if key in self.parameters:
                    savebin(model[key][iproc], path, iproc, prefix+key+suffix)
                else:
                    src = PATH.OUTPUT +'/'+ 'model_init'
                    dst = path
                    copybin(src, dst, iproc, prefix+key+suffix)

            if 'rho' in self.parameters:
                savebin(model['rho'][iproc], path, iproc, prefix+'rho'+suffix)
            else:
                src = PATH.OUTPUT +'/'+ 'model_init'
                dst = path
                copybin(src, dst, iproc, 'rho')

