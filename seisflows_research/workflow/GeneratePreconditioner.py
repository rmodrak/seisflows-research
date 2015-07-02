
from glob import glob

import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError

from seisflows.seistools import adjoint

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver
import preprocess
import postprocess


class GeneratePreconditioner(object):

    def check(self):
        """ Checks parameters and paths
        """
        # check paths
        if 'GLOBAL' not in PATH:
            raise ParameterError(PATH, 'GLOBAL')

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'OUTPUT' not in PATH:
            raise ParameterError(PATH, 'OUTPUT')

        # check input
        if 'DATA' not in PATH:
            setattr(PATH, 'DATA', None)

        if 'MODEL_INIT' not in PATH:
            raise ParameterError(PATH, 'MODEL_INIT')

        if 'CLIP' not in PAR:
            setattr(PAR, 'CLIP', 0.)

        # assertions
        0. <= PAR.CLIP <= 1.


    def main(self):
        # prepare directory structure
        unix.rm(PATH.GLOBAL)
        unix.mkdir(PATH.GLOBAL)

        # set up workflow machinery
        preprocess.setup()
        postprocess.setup()

        print 'Generating preconditioner...'
        system.run('solver', 'generate_preconditioner',
                   hosts='all',
                   process_traces=process_traces,
                   model_path=PATH.MODEL_INIT,
                   model_name='model',
                   model_type='gll')

        postprocess.combine_kernels(
            path=PATH.GLOBAL)

        self.process_kernels(
            path=PATH.GLOBAL)

        # save preconditioner
        src = PATH.GLOBAL +'/'+ 'kernels/absval'
        dst = PATH.OUTPUT +'/'+ 'precond'
        unix.cp(src, dst)

        print 'Finished\n'


    def process_kernels(self, path):
        assert (exists(path))

        # take absolute value
        parts = solver.load(path +'/'+ 'kernels/sum')
        for key in solver.parameters:
            parts[key] = np.abs(parts[key])

        solver.save(path +'/'+ 'kernels/absval',
                    parts,
                    suffix='_kernel')

        # clip
        if PAR.CLIP > 0.:
            for key in solver.parameters:
                maxval = np.max(parts[key])
                if maxval > 0:
                    parts[key] = np.clip(parts[key], 0., PAR.CLIP*maxval)

            src = PATH.GLOBAL +'/'+ 'kernels/absval'
            dst = PATH.GLOBAL +'/'+ 'kernels/absval_noclip'
            unix.mv(src, dst)

            solver.save(path +'/'+ 'kernels/absval',
                        parts,
                        suffix='_kernel')

        # smooth
        if PAR.SMOOTH > 0.:
            system.run('solver', 'smooth',
                       hosts='head',
                       path=path +'/'+ 'kernels/absval',
                       span=PAR.SMOOTH)

        # normalize
        parts = solver.load(path +'/'+ 'kernels/absval')
        for key in solver.parameters:
            parts[key] = np.mean(parts[key])/parts[key]
            solver.save(path +'/'+ 'kernels/absval',
                        parts,
                        suffix='_kernel')


def process_traces(path):
    unix.cd(path)

    s, h = preprocess.load(prefix='traces/syn/')
    s = preprocess.apply(preprocess.process_traces, [s], [h])

    s = preprocess.apply(adjoint.precond, [s], [h])
    preprocess.save(s, h, prefix='traces/adj/')

