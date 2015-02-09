import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import loadtxt, savetxt
from seisflows.tools.config import loadclass, ParameterObj
from seisflows.optimize import lib

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class Newton(loadclass('preprocess', 'base')):

    def prepare_apply_hess(self, path='.'):
        """ Prepares solver to compute action of Hessian
        """
        unix.cd(path)

        d, h = self.load(prefix='traces/obs/')
        s, _ = self.load(prefix='traces/lcg/')

        d = self.apply(self.process_traces, [d], [h])
        s = self.apply(self.process_traces, [s], [h])

        self.apply(self.write_residuals, [s, d], [h], inplace=False)

        s = self.apply(self.generate_adjoint_traces, [s, d], [h])
        self.save(s, h, prefix='traces/adj/')

