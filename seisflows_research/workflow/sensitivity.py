
from seisflows.tools.config import SeisflowsParameters, SeisflowsPaths, \
    ParameterError, loadclass

PAR = SeisflowsParameters()
PATH = SeisflowsPaths()

import system
import solver


class SensitivityAnalysis(object):
    """ Tests solver by running forward simulation
    """

    def check(self):
        """ Checks parameters and paths
        """

        # check paths
        if 'GLOBAL' not in PATH:
            raise Exception

        if 'LOCAL' not in PATH:
            setattr(PATH, 'LOCAL', None)

        if 'MODEL_INIT' not in PATH:
            raise Exception

        if 'MODEL_PERT' not in PATH:
            raise Exception

        if 'OUTPUT' not in PATH:
            raise Exception


    def main(self):
        """ Generates seismic data
        """

        print 'Running solver...'

        system.run('solver', 'generate_data',
                   hosts='all',
                   model_path=PATH.MODEL,
                   model_type='gll',
                   model_name='model')

        print "Finished"
