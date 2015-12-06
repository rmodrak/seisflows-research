
import numpy as np

from seisflows.tools import unix
from seisflows.tools.code import Struct, exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

from seisflows.seistools import adjoint, misfit, sbandpass, smute, readers, writers

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


class weighted(loadclass('preprocess', 'base')):
    """ Data preprocessing class
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(weighted, self).check()

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
        super(weighted, self).setup()


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

        # load source weights
        if PATH.SOURCE_WEIGHTS:
            ws = np.loadtxt(PATH.SOURCE_WEIGHTS)
        else:
            ws = np.ones(PAR.NSRC)

        # load receiver weights
        if PATH.RECEIVER_WEIGHTS:
            wr = np.loadtxt(PATH.RECEIVER_WEIGHTS)
        else:
            wr = np.ones(PAR.NREC)

        if PAR.EXPONENT:
            wr **= PAR.EXPONENT
            ws **= PAR.EXPONENT

        # apply source weights
        for key in self.channels:
            s[key] *= ws[system.getnode()]
            r[key] *= ws[system.getnode()]

        # apply receiver weights
        for key in self.channels:
            for ir in range(h.nr):
                if wr.ndim == 1:
                    s[key][:,ir] *= wr[ir]
                    r[key][ir] *= wr[ir]
                elif wr.ndim == 2:
                    s[key][:,ir] = np.dot(s[key][:,:],wr[:,ir])
                    r[key][ir] *= wr[ir,ir]

        self.save(s, h, prefix='traces/adj/')

        for key in self.channels:
            np.savetxt('residuals', r[key])


