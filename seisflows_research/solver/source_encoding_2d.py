
import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, custom_import

from seisflows.seistools.shared import SeisStruct
import seisflows.seistools.specfem2d as solvertools

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import preprocess


class source_encoding_2d(custom_import('solver', 'specfem2d_legacy')):
    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(source_encoding_2d, self).check()


    def initialize_solver_directories(self):
        """ Sets up directory in which to run solver
        """
        if 'NT_PADDED' not in PAR:
            raise Exception

        super(source_encoding_2d, self).initialize_solver_directories()
        solvertools.setpar('NSOURCES', PAR.NSRC)
        solvertools.setpar('nt', PAR.NT_PADDED)


    def write_receivers(self):
        """ Writes receivers file
        """
        unix.cd(self.getpath)

        _, hdr = preprocess.load('traces/obs')
        solvertools.write_receivers(hdr)

    def write_sources(self, sinfo=[], mapping=lambda i: [i]):
        """ Writes sources file
        """
        unix.cd(self.getpath)

        nodes = mapping(system.getnode())
        lines = []
        for i in nodes:
            solvertools.write_sources(vars(PAR), sinfo[i])
            with open('DATA/SOURCE', 'r') as f:
                lines.extend(f.readlines())

        with open('DATA/SOURCE', 'w') as f:
            f.writelines(lines)

    def initialize_adjoint_traces(self):
        zeros = np.zeros((PAR.NT_PADDED, PAR.NREC))
        h = SeisStruct(PAR.NREC, PAR.NT_PADDED, PAR.DT)
        for channel in ['x', 'y', 'z']:
            preprocess.writer(zeros, h, channel=channel, prefix='traces/adj')

