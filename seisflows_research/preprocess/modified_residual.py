
import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import Struct, exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

from seisflows.seistools import adjoint, misfit, sbandpass, smute, readers, writers

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class modified_residual(loadclass('preprocess', 'base')):
    """ Data preprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(modified_residual, self).check()

        if 'EXPONENT' not in PAR:
            setattr(PAR, 'EXPONENT', 1.0)

        if 'RECEIVER_WEIGHTS' not in PATH:
            setattr(PATH, 'RECEIVER_WEIGHTS', None)

        if 'SOURCE_WEIGHTS' not in PATH:
            setattr(PATH, 'SOURCE_WEIGHTS', None)

        if PATH.RECEIVER_WEIGHTS:
            assert exists(PATH.RECEIVER_WEIGHTS)

        if PATH.SOURCE_WEIGHTS:
            assert exists(PATH.SOURCE_WEIGHTS)


    def setup(self):
        """ Called once at beginning of workflow to perform any required
          initialization or setup tasks
        """
        super(modified_residual, self).setup()


    def load_weights(self, h):
        if PATH.RECEIVER_WEIGHTS:
            wr = np.loadtxt(PATH.RECEIVER_WEIGHTS)
        else:
            wr = self.compute_weights(h.rx, h.rz)

        if PATH.SOURCE_WEIGHTS:
            ws = np.loadtxt(PATH.SOURCE_WEIGHTS)
        else:
            ws = self.compute_weights(h.sx, h.sz)

        if PAR.EXPONENT:
            wr **= PAR.EXPONENT
            ws **= PAR.EXPONENT

        self.weights_rec = wr
        self.weights_src = ws


    def prepare_eval_grad(self, path='.'):
        """ Prepares solver for gradient evaluation by writing residuals and
          adjoint traces
        """
        import system

        unix.cd(path)

        d, h = self.load(prefix='traces/obs/')
        s, _ = self.load(prefix='traces/syn/')

        d = self.apply(self.process_traces, [d], [h])
        s = self.apply(self.process_traces, [s], [h])

        r = self.apply(self.write_residuals, [s, d], [h], inplace=False)
        s = self.apply(self.generate_adjoint_traces, [s, d], [h])

        if not hasattr(self, 'weights_rec') or \
           not hasattr(self, 'weights_src'): self.load_weights(h)

        wr = self.weights_rec
        ws = self.weights_src

        # apply receiver weights
        for key in self.channels:
            for ir in range(h.nr):
                if wr.ndim == 1:
                    s[key][:,ir] *= wr[ir]
                    r[key][ir] *= wr[ir]
                elif wr.ndim == 2:
                    s[key][:,ir] = np.dot(s[key][:,:],wr[:,ir])
                    r[key][ir] *= wr[ir,ir]

        # apply source weights
        for key in self.channels:
            s[key] *= ws[system.getnode()]
            r[key] *= ws[system.getnode()]

        self.save(s, h, prefix='traces/adj/')

        for key in self.channels:
            np.savetxt('residuals', r[key])


    def compute_weights(self, x, y):
        L = PAR.LENGTH
        n = x.size
        w = np.zeros(n)
        for i in range(n):
            R = distance(x[i], y[i], x, y)
            w[i] = sum(np.exp(-(R/L)**2.))**(-1.)
        w /= sum(w)
        w *= n
        return w


def distance(lat1, lon1, lat2, lon2):
    """ Haversine formula, modified from script by Wayne Dyck
    """
    dlat = np.radians(lat2-lat1)
    dlon = np.radians(lon2-lon1)
    a = np.sin(dlat/2) * np.sin(dlat/2) + np.cos(np.radians(lat1)) \
        * np.cos(np.radians(lat2)) * np.sin(dlon/2) * np.sin(dlon/2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    return c

