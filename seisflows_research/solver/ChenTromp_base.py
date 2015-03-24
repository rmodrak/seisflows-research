from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class ChenTromp_base(loadclass('solver', 'ChenTromp_base')):

    # model parameters expected by solver
    solver_parameters = []
    solver_parameters += ['rho']
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


    def export_kernels(self, path):
        super(specfem3d_ChenTromp, self).export_kernels(path)
        try:
            name = 'azimuth'
            src = join(glob(self.databases +'/'+ '*'+ name+'.bin'))
            dst = join(path, 'azimuth')
            unix.mv(src, dst)
        except:
            pass

