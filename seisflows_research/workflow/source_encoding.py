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


class source_encoding(custom_import('workflow', 'inversion')):
    """ Source encoding subclass
    """
    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(source_encoding, self).check()

        # check source encoding parameters
        if 'ENCODING' not in PAR:
            setattr(PAR, 'ENCODING', 1)

        if 'FAILRATE' not in PAR:
            setattr(PAR, 'FAILRATE', 0.)

        if 'SHIFT' not in PAR and PAR.ENCODING in [3,4]:
            raise Exception

        if 'NT_PADDED' not in PAR:
            if PAR.ENCODING in [3,4]:
                PAR.NT_PADDED = PAR.NT + (PAR.NSRC-1)*PAR.SHIFT
            else:
                PAR.NT_PADDED = PAR.NT


    def setup(self):
        """ Lays groundwork for inversion
        """
        super(source_encoding, self).setup()
        self.prepare_data()


    def initialize(self):
        """ Prepares for next model update iteration
        """
        stats = {}
        stats['ws'] = self.prepare_sources()
        stats['wr'] = self.prepare_receivers()

        # combine observations into 'supergather'
        self.combine(stats, tag='obs')

        # update input files
        solver.write_receivers()

        solver.write_sources(
            stats=stats,
            mapping=lambda _: range(PAR.NSRC))

        # generate synthetics
        super(source_encoding, self).initialize()


    def prepare_sources(self):
        """ Generates source encoding factors
        """
        if PAR.ENCODING == 0:
            ts = np.zeros(PAR.NSRC)
            fs = np.ones(PAR.NSRC)

        elif PAR.ENCODING == 1:
            # binary weights
            ts = np.zeros(PAR.NSRC)
            fs = np.sign(np.random.randn(PAR.NSRC))

        elif PAR.ENCODING == 2:
            # Gaussian weights
            ts = np.zeros(PAR.NSRC)
            fs = np.random.randn(PAR.NSRC)

        elif PAR.ENCODING == 3:
            # "staggered" shifts
            ts = np.arange(0, PAR.NSRC)*PAR.SHIFT*PAR.DT
            random.shuffle(ts)
            fs = np.ones(PAR.NSRC)

        elif PAR.ENCODING == 4:
            # random shifts
            ts = ((PAR.NSRC - 1)*PAR.SHIFT*PAR.DT)*np.random.rand(PAR.NSRC)
            fs = np.ones(PAR.NSRC)

        # collect factors
        stats = {}
        for i in range(PAR.NSRC):
            s[i].ts = ts[i]
            s[i].fs = fs[i]
        return stats


    def prepare_receivers(self):
        """ Generates receiver factors
        """
        if optimize.iter == 1:
            # generate receiver factors
            if PAR.FAILRATE == 0:
                rs = np.ones((PAR.NREC, PAR.NSRC))
            else:
                rs = np.random.rand(PAR.NREC, PAR.NSRC)
                rs = (rs > PAR.FAILRATE).astype(int)
            np.savetxt(PATH.GLOBAL + '/' + 'rs', rs, '%3d')

        # collect factors
        r = []
        rs = np.loadtxt(PATH.GLOBAL + '/' + 'rs')
        for i in range(PAR.NSRC):
            r.append(rs[:, i])
        return r


    def combine(self, stats, tag='obs'):
        """ Combines multiple sources into "supergather"
        """
        data = self.prepare_data(tag)

        for filename in self.filenames:
            summed = Stream()

            for dirname in self.dirnames:
                stream = copy(data[dirname][filename])

                # calculate time offset
                imin = int(stats.ts/PAR.DT)
                imax = imin + PAR.NT

                for ir in range(PAR.NREC):
                    # apply weights
                    stream[ir] *= sinfo[ir]
                    stream[ir] *= rinfo[ir]

                    summed[ir].data[imin:imax] += stream[ir].data[imin:imax]

            # save results
            preprocess.writer('', filename, summed)


    def prepare_data(self, tag='obs'):
        """ Loads from multiple sources into memory
        """
        # check if data already in memory
        if tag in globals():
            return globals()[tag]

        for dirname in self.dirnames:
            fullpath = PATH.DATA +'/'+ ''
            for filename in self.filenames:
                data[dirname][filename] = preprocess.reader(fullpath, filename)
        globals()[tag] = data

        return data


    @property
    def dirnames(self):
        return solver.check_source_names[0:PAR.NSRC]

    @property
    def filenames(self):
        return solver.data_filenames



def cdiff_adjoint(wsyn, wobs, nt, dt):
    # cross correlation difference
    cdiff = _np.correlate(wobs,wsyn) - _np.correlate(wobs,wobs)
    wadj = _np.convolve(wobs,cdiff)
    return 1e-10 * wadj


def cdiff_misfit(wsyn, wobs, nt, dt):
    cdiff = np.correlate(wobs, wsyn) - np.correlate(wobs, wobs)
    return np.sqrt(np.sum(cdiff*cdiff*dt))


