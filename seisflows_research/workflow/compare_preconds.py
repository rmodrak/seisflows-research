
from glob import glob

import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import cast, exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

from seisflows.seistools import adjoint
from seisflows.seistools.io import loadbypar, loadbin, savebin


PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver
import preprocess
import postprocess

migration = loadclass('workflow','migration')()



class compare_preconds(loadclass('workflow', 'compute_precond')):

    def check(self):
        """ Checks parameters and paths
        """
        migration.check()        
        super(compare_preconds, self).check()

        if 'CLIP' not in PAR:
            setattr(PAR, 'CLIP', True)


    def main(self):
        path = PATH.GLOBAL

        unix.rm(path)
        unix.mkdir(path)

        # set up workflow machinery
        preprocess.setup()
        postprocess.setup()

        print 'Generating preconditioner...'
        system.run('solver', 'generate_precond',
                   hosts='all',
                   process_traces=process_traces,
                   model_path=PATH.MODEL_INIT,
                   model_name='model',
                   model_type='gll')

        postprocess.combine_kernels(path,
            solver.parameters)

        for span in cast(PAR.SMOOTH):
            self.process_kernels(path, solver.parameters, 'precond',
                span=span)

            # save preconditioners
            names = glob(path +'/'+ 'precond')
            for src in names:
                dst = (PATH.OUTPUT +'/'+ '%s_%04d') % (unix.basename(src), span)
                unix.mv(src, dst)

        self.combine_kernels_hessian(path)
        self.process_kernels_hessian(path)

        print 'Finished\n'


    def combine_kernels_hessian(self, path):
        postprocess.combine_kernels(path,
             ['hessian1'] +
             ['hessian2'])

        fullpath = path +'/'+ 'kernels/sum'
        for iproc in [0]:#range(solver.mesh.nproc):
           hessian1 = loadbin(fullpath, iproc, 'hessian1_kernel')
           hessian2 = loadbin(fullpath, iproc, 'hessian2_kernel')
           savebin(hessian1 + hessian2, fullpath, iproc, 'hessian3_kernel')


    def process_kernels_hessian(self, path):
        for span in cast(PAR.SMOOTH):
            self.process_kernels(path, ['hessian1'], 'hessian1',
                span=span)

            self.process_kernels(path, ['hessian2'], 'hessian2',
                span=span)

            self.process_kernels(path, ['hessian3'], 'hessian3',
                span=span)

            # save preconditioners
            names = glob(path +'/'+ 'hessian?')
            for src in names:
                dst = (PATH.OUTPUT +'/'+ '%s_%04d') % (unix.basename(src), span)
                unix.mv(src, dst)


    def process_kernels(self, path, parameters, tag, span=0.):
        assert exists(path)
        assert len(parameters) > 0

        parts = self.load(
            path +'/'+ 'kernels/sum', 
            parameters)

        # clip
        if PAR.CLIP:
            for key in parameters:
                #parts[key] = np.abs(parts[key])
                parts[key][0] = np.clip(parts[key][0], 0, np.inf)
            self.save(path +'/'+ tag+'_nosmooth', parts)

        # smooth
        for key in parameters:
            parts[key][0] = self.smooth(parts[key][0], span)

        # normalize
        for key in parameters:
            meanval = np.mean(parts[key][0])
            parts[key][0] /= meanval
            #parts[key][0] **= -1
        self.save(path +'/'+ tag, parts)


    def save(self, path, parts, prefix='', suffix=''):
        unix.rm(path)
        unix.cp(PATH.MODEL_INIT, path)

        keys = parts.keys()
        assert len(keys) == 1

        precond_key = parts.keys()[0]
        for iproc in range(solver.mesh.nproc):
            for model_key in solver.parameters:
                savebin(parts[precond_key][iproc], 
                        path, iproc, prefix+model_key+suffix)


    def load(self, path, parameters, prefix='', suffix='_kernel'):
        parts = {}
        for key in parameters:
            parts[key] = []

        for iproc in range(solver.mesh.nproc):
            # read database files
            keys, vals = loadbypar(path, parameters, iproc, prefix, suffix)
            for key, val in zip(keys, vals):
                parts[key] += [val]
        return parts


    def smooth(self, v, span=0.):
        """ Smooths SPECFEM2D kernels by convolving them with a Gaussian
        """
        from seisflows.tools.array import meshsmooth, stack
        assert solver.mesh.nproc == 1

        if not span:
            return v

        # set up grid
        _,x = loadbypar(PATH.MODEL_INIT, ['x'], 0)
        _,z = loadbypar(PATH.MODEL_INIT, ['z'], 0)
        mesh = stack(x[0], z[0])

        vsmooth = meshsmooth(v.T, mesh, span)
        return vsmooth



def process_traces(path):
    unix.cd(path)

    s, h = preprocess.load(prefix='traces/syn/')

    s = preprocess.apply(adjoint.precond, [s], [h])

    preprocess.save(s, h, prefix='traces/adj/')


