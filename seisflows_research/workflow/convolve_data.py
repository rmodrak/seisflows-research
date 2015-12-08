
import numpy as np

from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

from seisflows.tools import unix
from seisflows.seistools import wavelets
from seisflows.seistools.signal import sbandpass, sconvolve

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver
import preprocess


class convolve_data(loadclass('workflow', 'test_forward')):
    """ Computes Green's functions, then convolves them with source time
      function after the fact
    """

    def check(self):
        """ Checks parameters and paths
        """
        super(generate_data, self).check()

        assert 'DT' in PAR
        assert 'NT' in PAR
        assert 'F0' in PAR


    def main(self):
        unix.rm(PATH.GLOBAL)
        unix.mkdir(PATH.GLOBAL)

        preprocess.setup()

        super(generate_data, self).main()

        print 'Convolving wavelet...'

        system.run('workflow', 'convolve_wavelet',
            hosts='all')

        print 'Finished\n'


    def convolve_wavelet(self):

        path = solver.getpath
        unix.cd(path)

        # load data
        s,h = preprocess.load('traces/obs')

        # convolve wavelet 
        if PAR.WAVELET:
            w = self.wavelet()
            s = preprocess.apply(sconvolve, [s], [h, w])

        if PAR.FREQLO or PAR.FREQHI:
            s = preprocess.apply(sbandpass, [s], [h, PAR.FREQLO, PAR.FREQHI])

        # save data
        dirname = 'traces_convolved'
        unix.mkdir(dirname)
        preprocess.save(s, h, dirname, suffix='.su')


    def wavelet(self):
        return getattr(wavelets, PAR.WAVELET)(PAR.DT, PAR.FP)
