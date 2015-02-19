
from seisflows.seistools.io import savebin, applymap

from seisflows.tools import unix
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')


class specfem3d_elastic(loadclass('solver', 'specfem3d')):
    """ Python interface for SPECFEM2D

      See base class for method descriptions
    """

    if PAR.MATERIALS == 'rho_bulk_c_bulk_beta':
        model_parameters = []
        model_parameters += ['vp']
        model_parameters += ['vs']

        kernel_parameters = []
        kernel_parameters += ['bulk_c']
        kernel_parameters += ['bulk_beta']

        raise NotImplementedError


    elif PAR.MATERIALS == 'rho_kappa_mu':
        model_parameters = []
        model_parameters += ['rho']
        model_parameters += ['vp']
        model_parameters += ['vs']

        model_parameter_maps = []
        model_parameter_maps += [['rho', lambda rho,kappa,mu : rho]]
        model_parameter_maps += [['vp', lambda rho,kappa,mu : ((kappa+1.333*mu)/rho)**0.5]]
        model_parameter_maps += [['vs', lambda rho,kappa,mu : (mu/rho)**0.5]]

        kernel_parameters = []
        kernel_parameters += ['rho']
        kernel_parameters += ['kappa']
        kernel_parameters += ['mu']

        kernel_parameter_maps = []
        kernel_parameter_maps += [['rho', lambda rho,vp,vs : rho]]
        kernel_parameter_maps += [['kappa', lambda rho,vp,vs : rho*(vp**2-1.333*vs**2)]]
        kernel_parameter_maps += [['mu', lambda rho,vp,vs : rho*vs**2]]


    def load(self, path, suffix='', verbose=False, **dummy):
        """ reads SPECFEM model
        """
        parameters = self.model_parameters
        maps = self.kernel_parameter_maps

        model = super(specfem3d_elastic, self).load(path, parameters, maps, suffix, verbose)
        return model


    def save(self, path, model):
        unix.mkdir(path)

        for iproc in range(PAR.NPROC):
            keys = []
            vals = []
            for key in self.kernel_parameters:
                keys += [key]
                vals += [model[key][iproc]]

            for key,val in zip(*applymap(vals, self.model_parameter_maps)):
                savebin(val, path, iproc, key)
            
    @property
    def parameters(self):
        return self.kernel_parameters
