
from glob import glob
from os.path import join

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


class compute_sensitivity(object):

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

        if 'MODEL' not in PATH:
            raise ParameterError(PATH, 'MODEL')

        if 'PERTURB' not in PATH:
            raise ParameterError(PATH, 'PERTURB')

        if 'MODEL_INIT' not in PATH:
            setattr(PATH, 'MODEL_INIT', PATH.MODEL)

        if 'MODEL_TRUE' not in PATH:
            setattr(PATH, 'MODEL_TRUE', PATH.MODEL)


    def main(self):
        """ Computes point spread from a peturbation dm with respect to a
          refrence model m
        """
        self.setup()

        m = solver.merge(solver.load(PATH.MODEL))
        dm = solver.merge(solver.load(PATH.PERTURB))

        self.evaluate_gradient(m+dm, join(PATH.GLOBAL, 'eval1'))
        self.evaluate_gradient(m-dm, join(PATH.GLOBAL, 'eval2'))


    def setup(self):
        path = PATH.GLOBAL

        # prepare directory structure
        unix.rm(path)
        unix.mkdir(path)

        # set up workflow machinery
        preprocess.setup()
        postprocess.setup()

        system.run('solver', 'setup',
                   hosts='all')


    def evaluate_gradient(self, model, path):
        """ Performs forward simulation to evaluate objective function
        """
        solver.save(join(path, 'model'), solver.split(model))

        system.run('solver', 'eval_func',
                   hosts='all',
                   path=path)

        system.run('solver', 'eval_grad',
                   hosts='all',
                   path=path)

        postprocess.write_gradient(
            path=path)

