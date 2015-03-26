
from glob import glob
from os.path import join

from seisflows.tools import unix

from seisflows.seistools.io import copybin, savebin
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system


class ChenTromp_base(loadclass('solver', 'specfem3d_legacy')):

    # model parameters expected by solver
    solver_parameters = []
    solver_parameters += ['A']
    solver_parameters += ['C']
    solver_parameters += ['N']
    solver_parameters += ['L']
    solver_parameters += ['F']
    solver_parameters += ['Jc']
    solver_parameters += ['Js']
    solver_parameters += ['Kc']
    solver_parameters += ['Ks']
    solver_parameters += ['Mc']
    solver_parameters += ['Ms']
    solver_parameters += ['Gc']
    solver_parameters += ['Gs']
    solver_parameters += ['Bc']
    solver_parameters += ['Bs']
    solver_parameters += ['Hc']
    solver_parameters += ['Hs']
    solver_parameters += ['Dc']
    solver_parameters += ['Ds']
    solver_parameters += ['Ec']
    solver_parameters += ['Es']


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


    def export_kernels(self, path):
        super(ChenTromp_base, self).export_kernels(path)
        try:
            name = 'azimuth'
            src = join(glob(self.databases +'/'+ '*'+ name+'.bin'))
            dst = join(path, 'azimuth')
            unix.mv(src, dst)
        except:
            pass


    def export_model(self, path):
        if system.getnode() == 0:
            parameters = self.solver_parameters
            parameters += ['rho']

            for parameter in parameters:
                unix.mkdir(path)
                src = glob(join(self.model_databases, '*'+parameter+'.bin'))
                dst = path
                unix.cp(src, dst)

