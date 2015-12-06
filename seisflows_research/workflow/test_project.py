
from os.path import join
import sys
import numpy as np

from seisflows.tools import msg
from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import divides, exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver
import optimize
import preprocess
import postprocess


class test_project(loadclass('workflow', 'inversion')):
    def check(self):
        """ Checks parameters and paths
        """
        super(test_project, self).check()

        if 'REFERENCE' not in PAR:
            setattr(PAR, 'REFERENCE', False)


    def write_model(self, path='', suffix=''):
        """ Writes model in format used by solver
        """
        unix.mkdir(path)
        src = PATH.OPTIMIZE +'/'+ 'm_' + suffix
        dst = path +'/'+ 'model'

        if PAR.REFERENCE:
            model = (
                solver.merge(solver.load(PATH.MODEL_INIT)) +
                postprocess.project_to_gll(
                    loadnpy(src)))

            solver.save(dst, solver.split(model))

        else:
            solver.save(dst,
                solver.split(
                    postprocess.project_to_gll(
                        loadnpy(src))))


    def save_gradient(self):
        src = join(PATH.GRAD, 'gradient')
        dst = join(PATH.OUTPUT, 'gradient_%04d' % optimize.iter)
        unix.mv(src, dst)


    def save_model(self):
        src = PATH.OPTIMIZE +'/'+ 'm_new'
        dst = join(PATH.OUTPUT, 'model_%04d' % optimize.iter)

        if PAR.REFERENCE:
            model = (
                solver.merge(solver.load(PATH.MODEL_INIT)) +
                postprocess.project_to_gll(
                    loadnpy(src)))

            solver.save(dst, solver.split(model))

        else:
            solver.save(dst,
                solver.split(
                    postprocess.project_to_gll(
                        loadnpy(src))))


    def save_kernels(self):
        src = join(PATH.GRAD, 'kernels')
        dst = join(PATH.OUTPUT, 'kernels_%04d' % optimize.iter)
        unix.mv(src, dst)

