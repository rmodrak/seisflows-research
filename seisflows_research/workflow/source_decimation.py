import random

import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import Struct, exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver
import optimize
import preprocess
import postprocess

from seisflows.workflow.inversion import base

class source_decimation(base):
    """ Source subset subclass
    """
    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(source_decimation, self).check()

        if 'NSRC_SUBSET' not in PAR:
            raise ParameterError

        assert PAR.NTASK == PAR.NSRC_SUBSET


    def setup(self):
        """ Lays groundwork for inversion
        """
        if PAR.BEGIN == 1:
            unix.rm(PATH.GLOBAL)
            unix.mkdir(PATH.GLOBAL)

            preprocess.setup()
            postprocess.setup()
            optimize.setup()

        for isrc in range(PAR.NSRC):
            solver.setup(subset=[isrc])


    def initialize(self):

        # choose subset
        solver.generate_subset()

        self.write_model(path=PATH.GRAD, suffix='new')

        print 'Generating synthetics'
        system.run('solver', 'eval_func',
                   hosts='all',
                   path=PATH.GRAD)

        self.sum_residuals(path=PATH.GRAD, suffix='new')

