
from os.path import basename, join

from seisflows.seistools.io import loadbypar

from seisflows.tools import unix
from seisflows.tools.code import Struct, exists
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()


def getstruct(args):
    return Struct(zip(*args))


def map(model, kernels):
    output = Struct()

    vp = model.vp
    vs = model.vs
    rho = model.rho
    kappa = rho*vp**2.
    mu = rho*vs**2.

    output.lame1 = [(1. - (2./3.)*(mu/kappa))*kernels.kappa]
    output.lame2 = [kernels.mu + (2./3.)*(mu/kappa)*kernels.kappa]
    output.rho = [rho]

    return output



class lambda_mu_2d(loadclass('solver', 'elastic2d')):
    """ Adds Lame parameter machinery to SPECFEM2D
    """
    assert PAR.MATERIALS == 'lambda_mu'

    def export_kernels(self, path):
        assert self.mesh.nproc == 1
        iproc = 0

        path = join(path, 'kernels')

        model_parameters = ['rho', 'vp', 'vs']
        kernel_parameters = ['rho', 'kappa', 'mu']

        model = getstruct(loadbypar(self.getpath+'/'+'DATA/', model_parameters, iproc))
        kernels = getstruct(loadbypar(self.getpath+'/'+'OUTPUT_FILES/', kernel_parameters, iproc, suffix='_kernel'))

        unix.mkdir(join(path, basename(self.path)))
        self.save(join(path, basename(self.path)), map(model, kernels), suffix='_kernel')

        

        
