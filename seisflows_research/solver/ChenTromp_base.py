
from glob import glob
from os.path import join

from seisflows.tools import unix

from seisflows.seistools.io import copybin, savebin
from seisflows.tools.config import loadclass, ParameterObj

PAR = ParameterObj('SeisflowsParameters')
PATH = ParameterObj('SeisflowsPaths')

import system


class ChenTromp_base(loadclass('solver', 'specfem3d_legacy')):

    # model parameters expected by solver
    solver_parameters = []
    solver_parameters += ['A']
    solver_parameters += ['C']
    solver_parameters += ['N']
    solver_parameters += ['L']
    solver_parameters += ['F']
    solver_parameters += ['Jc']
    solver_parameters += ['Js']
    solver_parameters += ['Kc']
    solver_parameters += ['Ks']
    solver_parameters += ['Mc']
    solver_parameters += ['Ms']
    solver_parameters += ['Gc']
    solver_parameters += ['Gs']
    solver_parameters += ['Bc']
    solver_parameters += ['Bs']
    solver_parameters += ['Hc']
    solver_parameters += ['Hs']
    solver_parameters += ['Dc']
    solver_parameters += ['Ds']
    solver_parameters += ['Ec']
    solver_parameters += ['Es']


    def save(self, path, model, prefix='', suffix=''):
        super(ChenTromp_base, self).save(path, model, prefix, suffix, solver_parameters)
