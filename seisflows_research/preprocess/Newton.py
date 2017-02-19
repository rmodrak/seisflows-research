
import sys
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.tools import loadtxt, savetxt
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class Newton(custom_import('preprocess', 'base')):

    def prepare_apply_hess(self, path='.'):
        """ Prepares solver to compute action of Hessian
        """
        solver = sys.modules['seisflows_solver']

        for filename in solver.data_filenames:
            obs = self.reader(path+'/'+'traces/obs', filename)
            syn = self.reader(path+'/'+'traces/lcg', filename)

            # process observations
            obs = self.apply_filter(obs)
            obs = self.apply_mute(obs)
            obs = self.apply_normalize(obs)

            # process synthetics
            syn = self.apply_filter(syn)
            syn = self.apply_mute(syn)
            syn = self.apply_normalize(syn)

            if PAR.MISFIT:
                self.write_residuals(path, syn, obs)

            self.write_adjoint_traces(path+'/'+'traces/adj', syn, obs, filename)


