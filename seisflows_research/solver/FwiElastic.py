
from os.path import join

from seisflows.seistools.io import copybin, loadbypar, savebin, splitvec, Minmax
from seisflows.seistools.io import Model as IOStruct

from seisflows.tools import unix
from seisflows.tools.code import exists
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')



class FwiElastic(object):
    """ Adds elastic inversion machinery
    """

    if 'MATERIALS' not in PAR:
        ParametersError(PAR, 'MATERIALS')

    if PAR.MATERIALS == 'bulk_c_bulk_mu':
        from seisflows.seistools.maps import bulk_c_bulk_mu_forward as map_forward
        from seisflows.seistools.maps import bulk_c_bulk_mu_inverse as map_inverse
        model_parameters = []
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['bulk_c']
        kernel_parameters += ['bulk_mu']

    elif PAR.MATERIALS == 'kappa_mu':
        from seisflows.seistools.maps import kappa_mu_forward as map_forward
        from seisflows.seistools.maps import kappa_mu_inverse as map_inverse
        model_parameters = []
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['kappa']
        kernel_parameters += ['mu']

    elif PAR.MATERIALS == 'vp_vs':
        from seisflows.seistools.maps import vp_vs_forward as map_forward
        from seisflows.seistools.maps import vp_vs_inverse as map_inverse
        model_parameters = []
        model_parameters += ['vp']
        model_parameters += ['vs']
        kernel_parameters = []
        kernel_parameters += ['vp']
        kernel_parameters += ['vs']

    if 'DENSITY' not in PAR:
        PAR.DENSITY = 'constant'

    if PAR.DENSITY == 'birch':
        import seisflows.seistools.maps.rho_birch as map_density

    elif PAR.DENSITY == 'constant':
        map_density = None

    elif PAR.DENSITY == 'variable':
        map_density = None
        model_parameters += ['rho']
        kernel_parameters += ['rho']



    def load(self, path, prefix='', suffix='', verbose=True):
        """ reads SPECFEM model or kernels
        """
        logpath = PATH.SUBMIT

        if 'kernel' in suffix:
            minmax = Minmax(self.kernel_parameters)

            # load kernels
            kernels = IOStruct(self.kernel_parameters)

            for iproc in range(PAR.NPROC):
                # read database files
                keys, vals = loadbypar(path, self.kernel_parameters, iproc, prefix, suffix)
                minmax.update(keys, vals)

                # store kernels
                for key, val in zip(keys, vals):
                    kernels[key] += [val]

            if verbose:
                minmax.write(path, logpath)
            return kernels

        else:
            minmax = Minmax(self.model_parameters)

            # load model
            model = IOStruct(self.kernel_parameters)

            for iproc in range(PAR.NPROC):
                # read database files
                keys, vals = loadbypar(path, self.model_parameters, iproc, prefix, suffix)
                minmax.update(keys, vals)

                if 'rho' not in keys:
                    key, val = loadbypar(path, ['rho'], iproc, prefix, suffix)
                    keys += key
                    vals += val

                # convert on the fly from one set of parameters to another
                mapped = self.map_forward(keys, vals)

                rho = mapped.pop('rho')
                if PAR.DENSITY == 'variable':
                     raise Exception
                for key, val in mapped.items():
                    model[key] += [val]

            if verbose:
                minmax.write(path, logpath)
            return model


    def save(self, path, obj, prefix='', suffix=''):
        unix.mkdir(path)
        model_init = join(PATH.OUTPUT, 'model_init')

        if 'kernel' in suffix:
            kernels = obj
            # write kernels
            for iproc in range(PAR.NPROC):
                keys = kernels.keys()
                vals = []
                for key in keys:
                    vals += [kernels[key][iproc]]

                # write database files
                for key, val in zip(keys, vals):
                    savebin(val, path, iproc, prefix+key+suffix)

        else:
            # write model
            model = obj
            for iproc in range(PAR.NPROC):
                keys = model.keys()
                vals = []
                for key in keys:
                    vals += [model[key][iproc]]
                if 'rho' not in keys:
                    key, val = loadbypar(model_init, ['rho'], iproc, prefix, suffix)
                    keys += key
                    vals += val

                # convert on the fly from one set of parameters to another
                mapped = self.map_inverse(keys, vals)

                # write database files
                rho = mapped.pop('rho')
                if PAR.DENSITY == 'variable':
                    savebin(rho, path, iproc, prefix+'rho'+suffix)
                elif PAR.DENSITY == 'constant':
                    savebin(rho, path, iproc, prefix+'rho'+suffix)
                else:
                    rho = self.map_density(keys, vals)
                    savebin(rho, path, iproc, prefix+'rho'+suffix)

                for key, val in mapped.items():
                    savebin(val, path, iproc, prefix+key+suffix)


    @property
    def parameters(self):
        return self.kernel_parameters
            

