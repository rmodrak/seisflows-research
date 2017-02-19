
import sys
import numpy as np

from os.path import exists
from obspy.core import Stream, Trace

from seisflows.plugins import adjoint, misfit
from seisflows.tools import unix
from seisflows.tools.tools import Struct
from seisflows.config import ParameterError, custom_import

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']

system = sys.modules['seisflows_system']


irec = 5

class double_difference(custom_import('preprocess', 'base')):
    """ Data preprocessing class
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        super(double_difference, self).check()

        if PAR.MISFIT not in ['Traveltime']:
            raise Exception


    def write_residuals(self, path, syn, dat):
        """ Computes residuals from observations and synthetics
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nr, _ = self.get_network_size(syn)

        d = np.zeros((nr, nr))
        dd = np.zeros((nr, nr))

        for i in range(nr):
            for j in range(irec):#range(i):

                print i,j
                d[i,j] = self.misfit_dd(syn[i].data,syn[j].data,nt,dt)
                dd[i,j] = d[i,j] - self.misfit_dd(dat[i].data,dat[j].data,nt,dt)

        # write residuals
        filename = path +'/'+ 'residuals'
        if exists(filename):
            rsd = list(np.loadtxt(filename))
        else:
            rsd = []
        for i in range(nr):
            rsd += [np.sum(abs(dd), 0)]
        np.savetxt(filename, rsd)

        np.savetxt(path +'/'+ 'dij', d)
        np.savetxt(path +'/'+ 'ddij', dd)


    def generate_adjoint_traces(self, path, syn, dat):
        """ Computes adjoint traces from observed and synthetic traces
        """
        nt, dt, _ = self.get_time_scheme(syn)
        nr, _ = self.get_network_size(syn)

        d = np.loadtxt(path +'/'+ 'dij')
        dd = np.loadtxt(path +'/'+ 'ddij')

        adj = Stream()
        for i in range(nr):
            adj.append(Trace(
                data=np.zeros(nt, dtype='float32'),
                header=syn[i].stats))

        for i in range(nr):
            for j in range(irec):#range(i):

                ai = adj[i].data
                aj = adj[j].data
                si = syn[i].data
                sj = syn[j].data

                tmpi = self.adjsrc_dd(sj, si, -d[i,j], nt, dt)
                tmpj = self.adjsrc_dd(si, sj, +d[i,j], nt, dt)

                ai += dd[i,j] * tmpi
                aj -= dd[i,j] * tmpj

        return adj



    def adjsrc_dd(self, si, sj, t0, nt, dt):
        w = np.zeros(nt)
        it = t0/dt
        print it

        w[1:-1] = (si[2:] - si[0:-2])/(2.*dt)

        if t0 > 0: w[it:] = w[:-it]
        if t0 < 0: w[-it:] = w[it:]

        w *= 1./(sum(w*sj*dt))


    def misfit_dd(self, si, sj, nt, dt):
        cc = abs(np.convolve(si, np.flipud(sj)))
        it = np.argmax(cc)
        return (it-nt+1)*dt

