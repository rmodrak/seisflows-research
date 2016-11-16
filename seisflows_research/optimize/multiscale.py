
from os.path import join
import sys
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.code import exists, loadtxt, savetxt
from seisflows.config import , \
    ParameterError, custom_import


PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

import solver


class multiscale(custom_import('optimize', 'debug')):
    """ Allows restarting between multiscale transitions
    """

    def check(self):
        super(multiscale, self).check()

        if not hasattr(PAR, 'NORESTART'):
            setattr(PAR, 'NORESTART', False)


    def compute_direction(self):
        """ Computes model update direction from stored gradient
        """
        self.multiscale_setup()
        if self.multiscale_status():
            self.multiscale_restart()
            return

        p_new = super(multiscale, self).compute_direction()

        if self.iter == PAR.BEGIN and PAR.NORESTART:
            self.restarted = True

        return p_new


    def multiscale_setup(self):
        if self.iter == PAR.BEGIN and self.iter > 1:
            if PAR.SCHEME in ['NLCG']:
                self.NLCG.path = PATH.OPTIMIZE
            elif PAR.SCHEME in ['LBFGS']:
                self.LBFGS.path = PATH.OPTIMIZE
            self.writer.path = PATH.OUTPUT
            self.stepwriter.path = PATH.SUBMIT


    def multiscale_status(self):
        if self.iter == 1:
            return False
        elif PAR.NORESTART:
            print 'MADE IT HERE'
            return False
        elif self.iter == PAR.BEGIN:
            return True
        else: 
            return False


    def multiscale_restart(self):
        self.multiscale_restarted = True
        self.stepwriter.iter += 1
        super(multiscale, self).restart()





