import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.tools import loadtxt, savetxt
from seisflows.config import , \
    ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class GaussNewton(custom_import('preprocess', 'base')):

    def prepare_apply_hess(self, path='.'):
        """ Prepares solver to compute action of Hessian
        """
        unix.cd(path)

        d, h = self.load(prefix='traces/syn/')
        s, _ = self.load(prefix='traces/lcg/')

        d = self.apply(self.process_traces, [d], [h])
        s = self.apply(self.process_traces, [s], [h])

        self.apply(self.write_residuals, [s, d], [h], inplace=False)

        s = self.apply(self.generate_adjoint_traces, [s, d], [h])
        self.save(s, h, prefix='traces/adj/')

